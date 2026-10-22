// swiftld - byte-identical, fast drop-in for vcftools 0.1.17 pairwise LD (--hap-r2 / --geno-r2).
//
// The specification is vcftools 0.1.17 itself (variant_file_output.cpp: output_haplotype_r2,
// output_genotype_r2, calc_hap_r2, calc_geno_r2; entry_filters.cpp; vcf_entry*.cpp; parameters.cpp).
// Every behaviour below that looks odd is there because the reference does it; comments cite the source.
//
// Licence: LGPL-3.0 (derived from vcftools semantics; see LICENSE_NOTE.txt).
#include <malloc.h>
#include <algorithm>
#include <atomic>
#include <cfloat>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <set>
#include <string>
#include <string_view>
#include <thread>
#include <chrono>
#include <functional>
#include <unordered_map>
#include <vector>

#include <fcntl.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "fmt.h"
#include "libdeflate/libdeflate.h"

using std::string;
using std::vector;
typedef std::string_view sv;

// ------------------------------------------------------------------------------------------------ errors
[[noreturn]] static void die(const string& msg, int code = 1) {
    fprintf(stderr, "Error: %s\n", msg.c_str());
    fflush(stderr);
    _exit(code);
}

// ------------------------------------------------------------------------------------------------ params
struct Params {
    string vcf_filename;
    bool vcf_format = false, compressed = false;
    string out = "out";
    bool stream_out = false;
    bool hap = false, geno = false;
    int num_outputs = 0;
    bool phased_only = false;
    int ld_bp_window_size = INT_MAX, ld_snp_window_size = INT_MAX;
    int ld_bp_window_min = -1, ld_snp_window_min = -1;
    double min_r2 = -1.0, min_maf = -1.0, max_maf = DBL_MAX, min_site_call_rate = 0;
    std::set<string> chrs;
    bool remove_filtered_all = false, remove_indels = false;
    int thin = -1;
    vector<string> keep_files, remove_files;
};

static Params parse_args(int argc, char** argv) {
    Params P;
    vector<string> a(argv, argv + argc);
    auto get_arg = [&](size_t i) -> const string& {
        if (i >= a.size()) die("Requested Missing Argument", 76);  // parameters::get_arg
        return a[i];
    };
    size_t i = 1;
    while (i < a.size()) {
        const string& s = a[i];
        if (s == "--vcf") { P.vcf_format = true; P.compressed = false; P.vcf_filename = get_arg(i + 1); i++; }
        else if (s == "--gzvcf") { P.vcf_format = true; P.compressed = true; P.vcf_filename = get_arg(i + 1); i++; }
        else if (s == "--out") { P.out = get_arg(i + 1); i++; }
        else if (s == "--stdout" || s == "-c") P.stream_out = true;
        else if (s == "--hap-r2") { P.hap = true; P.phased_only = true; P.num_outputs++; }
        else if (s == "--geno-r2") { P.geno = true; P.num_outputs++; }
        else if (s == "--ld-window-bp") { P.ld_bp_window_size = atoi(get_arg(i + 1).c_str()); i++; }
        else if (s == "--ld-window") { P.ld_snp_window_size = atoi(get_arg(i + 1).c_str()); i++; }
        else if (s == "--ld-window-bp-min") { P.ld_bp_window_min = atoi(get_arg(i + 1).c_str()); i++; }
        else if (s == "--ld-window-min") { P.ld_snp_window_min = atoi(get_arg(i + 1).c_str()); i++; }
        else if (s == "--min-r2") { P.min_r2 = atof(get_arg(i + 1).c_str()); i++; }
        else if (s == "--chr") { P.chrs.insert(get_arg(i + 1)); i++; }
        else if (s == "--maf") { P.min_maf = atof(get_arg(i + 1).c_str()); i++; }
        else if (s == "--max-missing") { P.min_site_call_rate = atof(get_arg(i + 1).c_str()); i++; }
        else if (s == "--remove-filtered-all") P.remove_filtered_all = true;
        else if (s == "--remove-indels") P.remove_indels = true;
        else if (s == "--thin") { P.thin = atoi(get_arg(i + 1).c_str()); i++; }
        else if (s == "--phased") P.phased_only = true;
        else if (s == "--keep") { P.keep_files.push_back(get_arg(i + 1)); i++; }
        else if (s == "--remove") { P.remove_files.push_back(get_arg(i + 1)); i++; }
        else if (s == "--temp") { get_arg(i + 1); i++; }  // no temporary files are used
        else die("Unsupported option for swiftld (LD subset only): " + s, 2);
        i++;
    }
    // parameters::check_parameters (the subset that can fire here)
    if (P.vcf_filename == "-") die("stdin input is not supported by swiftld", 2);
    if (P.num_outputs > 1) { fprintf(stderr, "Error: Only one output function may be called.\n"); exit(0); }
    if (P.vcf_filename == "") { fprintf(stderr, "Error: Input file required.\n"); exit(0); }
    if (P.min_maf != -1.0 && (P.min_maf < 0.0 || P.min_maf > 1.0)) die("MAF must be between 0 and 1.", 4);
    if (P.min_site_call_rate > 1) die("Minimum Call rate cannot be greater than 1.", 5);
    if (P.num_outputs == 0) die("swiftld needs --hap-r2 or --geno-r2", 2);
    return P;
}

// ------------------------------------------------------------------------------------------------ input
// Background decompression: the inflater thread fills the buffer front to back while the parser follows
// behind it. The vendored libdeflate is patched to publish, after every deflate block, how many output bytes
// are final (swiftld_progress_ptr); bytes below that offset are never written again.
extern "C" size_t* swiftld_progress_ptr;
extern "C" size_t swiftld_progress_base;
extern "C" void (*swiftld_block_hook)(size_t);
struct Stream {
    size_t progress = 0;   // written by the inflater with __atomic_store_n(release)
    std::atomic<bool> done{false}, failed{false};
    std::atomic<size_t> final_len{0};
    double t_inflate = 0;
    char* buf = nullptr;
    size_t cap = 0;
    // Page recycling: the parser publishes `safe`, the lowest offset it can still read (current line start, or
    // a retained field of an earlier line). After each deflate block the inflater moves the physical pages of
    // whole 2 MB chunks below `safe` (and > 128 KB behind its own output, beyond the 32 KB match window) to
    // 2 MB chunks it has not written yet, so most of the output lands in already-populated pages instead of
    // freshly faulted (zeroed, host-backed) ones. The old addresses become unmapped.
    std::atomic<size_t> safe{0};
    size_t recycle_next = 0;                  // inflater-owned: next chunk offset to recycle
    size_t recycled_until = 0;                // inflater-owned: targets below this offset are taken
    std::atomic<size_t> n_recycled{0};
    bool recycle_ok = true;
};
static Stream* g_recycle_stream = nullptr;
static void recycle_hook(size_t out_off) {
    Stream* st = g_recycle_stream;
    if (!st || !st->recycle_ok) return;
    const size_t CH = size_t(2) << 20;
    const size_t safe = st->safe.load(std::memory_order_acquire);
    const size_t behind = out_off > (size_t(128) << 10) ? out_off - (size_t(128) << 10) : 0;
    const size_t lim = std::min(safe, behind) & ~(CH - 1);
    size_t tgt = std::max(st->recycled_until, (out_off + (size_t(64) << 10) + CH - 1) & ~(CH - 1));
    while (st->recycle_next + CH <= lim && tgt + CH <= st->cap) {
        void* r = mremap(st->buf + st->recycle_next, CH, CH, MREMAP_MAYMOVE | MREMAP_FIXED, st->buf + tgt);
        if (r == MAP_FAILED) { st->recycle_ok = false; return; }
        st->recycle_next += CH;
        tgt += CH;
        st->recycled_until = tgt;
        st->n_recycled.fetch_add(1, std::memory_order_relaxed);
    }
}

struct Input {
    const char* data = nullptr;
    size_t len = 0;
    Stream* st = nullptr;   // non-null: data is being produced by a background inflater
    std::thread* st_thread = nullptr;
};
static std::atomic<bool> g_inflate_ok{false};          // background inflate finished without error
static std::atomic<bool> g_after_inflate_ready{false};
static std::function<void()> g_after_inflate_fn;        // what the inflater thread does next (pair work)

// Large anonymous buffer backed by transparent huge pages when the kernel allows (THP "madvise" mode):
// page-faulting 85 MB in 4 KB pages costs ~40 ms on this VM, ~15 ms in 2 MB pages.
static char* alloc_big(size_t n) {
    const size_t A = size_t(2) << 20;
    void* m = mmap(nullptr, n + A, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (m == MAP_FAILED) return (char*)malloc(n);
    char* a = (char*)(((uintptr_t)m + A - 1) & ~(uintptr_t)(A - 1));
    madvise(a, n, MADV_HUGEPAGE);
    return a;
}

// Bump allocator over one THP-advised reservation for long-lived buffers that are never freed (the output
// buffer pool, site blocks and their genotype planes). Faulting fresh memory costs ~1.4 us per 4 KB page on
// this VM against ~0.15-0.3 ms per 2 MB huge page (2.5-5x less per byte). Memory is zero (fresh anonymous
// pages); alloc returns nullptr when the reservation is exhausted and the caller falls back to calloc/malloc.
struct BumpArena {
    char* base = nullptr;
    size_t cap = 0;
    std::atomic<size_t> used{0};
    void init(size_t n) {
        void* m = mmap(nullptr, n + (size_t(2) << 20), PROT_READ | PROT_WRITE,
                       MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
        if (m == MAP_FAILED) return;
        base = (char*)(((uintptr_t)m + (size_t(2) << 20) - 1) & ~uintptr_t((size_t(2) << 20) - 1));
        madvise(base, n, MADV_HUGEPAGE);
        cap = n;
    }
    void* alloc(size_t n) {
        if (!base) return nullptr;
        n = (n + 63) & ~size_t(63);
        const size_t o = used.fetch_add(n, std::memory_order_relaxed);
        if (o + n > cap) return nullptr;
        return base + o;
    }
};
static BumpArena g_arena;

static int n_threads_env();

static Input read_input(const Params& P, bool allow_stream) {
    const string& fn = P.vcf_filename;
    struct stat st;
    if (stat(fn.c_str(), &st) != 0) die("Can't determine file type of " + fn);
    if (!S_ISREG(st.st_mode)) die("Does not appear to be a regular file: " + fn);
    if (fn.size() < 4) die("filename too short (vcftools throws std::out_of_range)");
    if (fn.substr(fn.size() - 4) == ".bcf") die("Filename ends in '.bcf'. Shouldn't you be using --bcf?");
    if (!P.compressed && fn.substr(fn.size() - 3) == ".gz")
        die("Filename ends in '.gz'. Shouldn't you be using --gzvcf or --gzdiff?");
    int fd = open(fn.c_str(), O_RDONLY);
    if (fd < 0) die("Could not open VCF file: " + fn);
    size_t flen = size_t(st.st_size);
    const unsigned char* raw = nullptr;
    if (flen) {
        void* m = mmap(nullptr, flen, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
        if (m == MAP_FAILED) die("mmap failed: " + fn);
        raw = (const unsigned char*)m;
    }
    close(fd);
    Input in;
    if (!P.compressed || flen < 2 || raw[0] != 0x1f || raw[1] != 0x8b) {
        // Plain text (zlib's gzread passes non-gzip data through unchanged).
        in.data = (const char*)raw;
        in.len = flen;
        return in;
    }
    // gzip: decompress all members (gzread semantics: concatenated members; trailing non-gzip ignored).
    // Fast path: libdeflate, one call per member into a buffer sized from the ISIZE trailer.
    size_t cap = 0;
    {
        uint32_t isize;
        memcpy(&isize, raw + flen - 4, 4);
        cap = size_t(isize) + 4096;
        if (cap < flen * 3) cap = flen * 3 + 65536;
    }
    char* out = alloc_big(cap);
    if (allow_stream) {
        // One attempt with libdeflate in a background thread. Any failure (output larger than the ISIZE-based
        // buffer, damaged member) sets `failed`; the parser then re-executes the program without streaming.
        Stream* st = new Stream();
        st->buf = out;
        st->cap = cap;
        g_recycle_stream = st;
        swiftld_block_hook = getenv("SWIFTLD_NORECYCLE") ? nullptr : recycle_hook;
        in.st_thread = new std::thread([st, raw, flen]() {
            auto t0 = std::chrono::steady_clock::now();
            libdeflate_decompressor* dec = libdeflate_alloc_decompressor();
            size_t olen = 0, pos = 0;
            bool ok = true;
            swiftld_progress_ptr = &st->progress;
            while (pos + 2 <= flen && raw[pos] == 0x1f && raw[pos + 1] == 0x8b) {
                size_t in_used = 0, out_got = 0;
                swiftld_progress_base = olen;
                libdeflate_result r = libdeflate_gzip_decompress_ex(dec, raw + pos, flen - pos, st->buf + olen,
                                                                    st->cap - olen, &in_used, &out_got);
                if (r != LIBDEFLATE_SUCCESS) { ok = false; break; }
                olen += out_got;
                pos += in_used;
            }
            libdeflate_free_decompressor(dec);
            swiftld_progress_ptr = nullptr;
            swiftld_block_hook = nullptr;
            st->final_len.store(olen, std::memory_order_relaxed);
            st->t_inflate = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            st->failed.store(!ok, std::memory_order_relaxed);
            if (ok) g_inflate_ok.store(true, std::memory_order_release);
            st->done.store(true, std::memory_order_release);
            if (!ok) return;  // the parser re-executes the program
            // become a pair worker once the main thread has set up the pair machinery
            while (!g_after_inflate_ready.load(std::memory_order_acquire)) sched_yield();
            g_after_inflate_fn();
        });
        in.data = out;
        in.len = 0;
        in.st = st;
        return in;
    }
    size_t olen = 0, pos = 0;
    bool ok = true;
    libdeflate_decompressor* dec = libdeflate_alloc_decompressor();
    while (pos + 2 <= flen && raw[pos] == 0x1f && raw[pos + 1] == 0x8b) {
        size_t in_used = 0, out_got = 0;
        libdeflate_result r;
        for (;;) {
            r = libdeflate_gzip_decompress_ex(dec, raw + pos, flen - pos, out + olen, cap - olen, &in_used, &out_got);
            if (r != LIBDEFLATE_INSUFFICIENT_SPACE) break;
            size_t ncap = cap * 2 + (1 << 20);
            char* nout = alloc_big(ncap);
            memcpy(nout, out, olen);
            out = nout;  // the old block is simply abandoned (process-lifetime memory)
            cap = ncap;
        }
        if (r != LIBDEFLATE_SUCCESS) { ok = false; break; }
        olen += out_got;
        pos += in_used;
    }
    libdeflate_free_decompressor(dec);
    if (!ok) {
        // corrupt or truncated member: redo with zlib streaming, which keeps the bytes before the damage
        // (what gzgets would have returned)
        olen = 0;
        pos = 0;
        while (pos + 2 <= flen && raw[pos] == 0x1f && raw[pos + 1] == 0x8b) {
            z_stream zs;
            memset(&zs, 0, sizeof zs);
            if (inflateInit2(&zs, 15 + 16) != Z_OK) die("zlib init failed");
            zs.next_in = (Bytef*)(raw + pos);
            zs.avail_in = uInt(std::min<size_t>(flen - pos, UINT_MAX));
            int ret;
            for (;;) {
                if (cap - olen < 65536) {
                    size_t ncap = cap * 2;
                    char* nout = alloc_big(ncap);
                    memcpy(nout, out, olen);
                    out = nout;
                    cap = ncap;
                }
                zs.next_out = (Bytef*)(out + olen);
                zs.avail_out = uInt(std::min<size_t>(cap - olen, UINT_MAX));
                uInt before = zs.avail_out;
                ret = inflate(&zs, Z_NO_FLUSH);
                olen += before - zs.avail_out;
                if (ret == Z_STREAM_END) break;
                if (ret == Z_DATA_ERROR || ret == Z_STREAM_ERROR || ret == Z_NEED_DICT || ret == Z_MEM_ERROR) {
                    // gzread reports the error and returns nothing more (gzgets -> NULL, gzeof stays false):
                    // vcftools then fails on the missing header, or loops forever reading empty lines.
                    die("compressed data error in gzip input (vcftools fails or does not terminate here)");
                }
                if (ret != Z_OK && ret != Z_BUF_ERROR) break;
                if (ret == Z_BUF_ERROR && zs.avail_in == 0) break;
            }
            pos += size_t(zs.next_in - (raw + pos));
            inflateEnd(&zs);
            if (ret != Z_STREAM_END) break;
        }
    }
    in.data = out;
    in.len = olen;
    return in;
}

// ------------------------------------------------------------------------------------------------ lines
// vcf_file::read_line: getline/gzgets, then erase trailing " \t\n\r".
static char** g_argv = nullptr;
[[noreturn]] static void restart_without_streaming() {
    setenv("SWIFTLD_NOSTREAM", "1", 1);
    execv("/proc/self/exe", g_argv);
    die("re-exec failed");
}

#ifndef MADV_POPULATE_WRITE
#define MADV_POPULATE_WRITE 23
#endif
struct LineReader {
    const char* p;
    const char* end;
    Stream* st = nullptr;       // background inflater, or null when the whole input is present
    std::function<bool()> idle; // work to do while waiting for decompressed data
    size_t prefault = 0;        // the buffer is populated ahead of the inflater up to this offset
    // While waiting and with no rows to format, fault in (allocate + zero) the output pages ahead of the
    // decompressor, so that the inflater thread writes into present pages: 85 MB decompress in 44 ms warm vs
    // 54-62 ms into fresh pages.
    bool no_prefault = false;   // set when the inflater recycles pages (prefaulted pages would be replaced)
    bool prefault_ahead() {
        if (no_prefault) return false;
        const size_t STEP = size_t(2) << 20, AHEAD = size_t(32) << 20;
        const size_t prog = __atomic_load_n(&st->progress, __ATOMIC_ACQUIRE);
        if (prefault < prog) prefault = (prog + STEP - 1) & ~(STEP - 1);
        if (prefault >= st->cap || prefault > prog + AHEAD) return false;
        const size_t n = std::min(STEP, st->cap - prefault);
        madvise(st->buf + prefault, n, MADV_POPULATE_WRITE);
        prefault += n;
        return true;
    }
    void refresh() {
        if (st->done.load(std::memory_order_acquire)) {
            if (st->failed.load(std::memory_order_relaxed)) restart_without_streaming();
            end = st->buf + st->final_len.load(std::memory_order_relaxed);
            st = nullptr;
            return;
        }
        const size_t prog = __atomic_load_n(&st->progress, __ATOMIC_ACQUIRE);
        if (st->buf + prog > end) end = st->buf + prog;
    }
    // Regular genotype lines without the newline search: find the 9th tab with 64-byte compares (stopping at a
    // newline) and accept the line end n4 bytes after it if that byte is a newline already in the buffer. The
    // end is then provisional (prov = true): the caller either proves the sample section has no newline (the
    // regular-line check accepts only genotype characters) or calls settle(). Otherwise this is next().
    bool next_pred(const char*& s, const char*& e, size_t n4, bool& prov) {
        prov = false;
#if defined(__AVX512F__) && defined(__AVX512BW__) && defined(__BMI2__)
        if (n4 && p < end) {
            const __m512i vt = _mm512_set1_epi8('\t'), vn = _mm512_set1_epi8('\n');
            const char* b = p;
            int tabs = 0;
            while (b < end && b - p < 4096) {
                const size_t rem = size_t(end - b);
                const __mmask64 lm = rem >= 64 ? ~__mmask64(0) : ((__mmask64(1) << rem) - 1);
                const __m512i v = _mm512_maskz_loadu_epi8(lm, b);   // bytes past `end` read as 0: no match
                const uint64_t mt = _mm512_cmpeq_epi8_mask(v, vt), mn = _mm512_cmpeq_epi8_mask(v, vn);
                const int c = __builtin_popcountll(mt);
                if (tabs + c >= 9) {
                    const uint64_t bit = _pdep_u64(uint64_t(1) << (8 - tabs), mt);   // the 9th tab of the line
                    if (mn & (bit - 1)) break;                                         // a newline comes first
                    const char* le = b + __builtin_ctzll(bit) + 1 + n4;
                    if (le < end && *le == '\n') { s = p; e = le; p = le + 1; prov = true; return true; }
                    break;
                }
                if (mn) break;
                tabs += c;
                b += 64;
            }
        }
#endif
        return next(s, e);
    }
    // The true end of a line returned provisionally as [s, e) (*e == '\n'): the first newline, then the same
    // trailing-whitespace trim as next(); the reader continues after that newline.
    void settle(const char* s, const char*& e) {
        const char* nl = (const char*)memchr(s, '\n', size_t(e - s) + 1);
        p = nl + 1;
        const char* le = nl;
        while (le > s && (le[-1] == ' ' || le[-1] == '\t' || le[-1] == '\n' || le[-1] == '\r')) le--;
        e = le;
    }
    bool next(const char*& s, const char*& e) {
        for (;;) {
            if (p < end) {
                const char* nl = (const char*)memchr(p, '\n', size_t(end - p));
                if (nl || !st) {
                    const char* le = nl ? nl : end;
                    s = p;
                    p = nl ? nl + 1 : end;
                    while (le > s && (le[-1] == ' ' || le[-1] == '\t' || le[-1] == '\n' || le[-1] == '\r')) le--;
                    e = le;
                    return true;
                }
            } else if (!st) {
                return false;
            }
            const char* old = end;
            refresh();
            if (st && end == old) {
                if (!(idle && idle()) && !prefault_ahead())
                    for (int k = 0; k < 64; k++) __builtin_ia32_pause();
            }
        }
    }
};

// split like `while(!ss.eof()) getline(ss, tok, delim)` (python str.split semantics)
static void split_sv(sv s, char d, vector<sv>& out) {
    out.clear();
    size_t st = 0;
    for (;;) {
        size_t k = s.find(d, st);
        if (k == sv::npos) { out.push_back(s.substr(st)); break; }
        out.push_back(s.substr(st, k - st));
        st = k + 1;
    }
}

// ------------------------------------------------------------------------------------------------ genotypes
// vcf_entry::set_indv_GENOTYPE_and_PHASE(indv, string) + header::str2int
static int str2int_allele(const char* s, size_t n) {
    if (n == 0) return -1;                      // header::str2int missing value
    if (n == 1 && s[0] == '.') return -1;       // caller tests in != "."
    char buf[64];
    size_t m = std::min<size_t>(n, 63);
    memcpy(buf, s, m);
    buf[m] = 0;
    return atoi(buf);
}

struct GT {
    int a1, a2;
    char phase;
    unsigned char ploidy;
    bool poly;
};

static inline void parse_gt(const char* g, size_t len, GT& o) {
    if (len == 3 && (g[1] == '/' || g[1] == '|')) {
        o.ploidy = 2;
        o.phase = g[1];
        o.a1 = (g[0] != '.') ? (g[0] - '0') : -1;
        o.a2 = (g[2] != '.') ? (g[2] - '0') : -1;
        o.poly = false;
        return;
    }
    size_t pos = len, pos2 = len;
    for (size_t k = 0; k < len; k++)
        if (g[k] == '/' || g[k] == '|') { if (pos == len) pos = k; pos2 = k; }
    if (pos != len) {
        o.ploidy = 2;
        o.phase = g[pos];
        sv f(g, pos), s(g + pos + 1, len - pos - 1);
        o.a1 = (f == ".") ? -1 : str2int_allele(f.data(), f.size());
        o.a2 = (s == ".") ? -1 : str2int_allele(s.data(), s.size());
    } else {
        o.ploidy = 1;
        o.phase = '|';
        sv f(g, len);
        o.a1 = (f == ".") ? -1 : str2int_allele(f.data(), f.size());
        o.a2 = -1;
    }
    o.poly = (pos != pos2);
}

// GT subfield of one sample column (vcf_entry::parse_genotype_entry with GT only)
static inline void parse_sample(const char* fs, const char* fe, int GT_idx, GT& o) {
    if (GT_idx < 0) { o.a1 = -1; o.a2 = -1; o.phase = '/'; o.ploidy = 2; o.poly = false; return; }
    int i = 0;
    const char* st = fs;
    for (;;) {
        const char* c = (const char*)memchr(st, ':', size_t(fe - st));
        if (!c) {
            if (i == GT_idx) { parse_gt(st, size_t(fe - st), o); return; }
            break;
        }
        if (i == GT_idx) { parse_gt(st, size_t(c - st), o); return; }
        i++;
        st = c + 1;
    }
    o.a1 = -1; o.a2 = -1; o.phase = '/'; o.ploidy = 2; o.poly = false;
}

// output_*_r2 temp-file byte encoding and read_temp_site decoding, reproduced exactly.
static inline void enc_dec(int first, int second, int& d1, int& d2) {
    char out_byte = 0x00;
    if (first == -1) out_byte |= 0x02; else out_byte |= (char)first;
    out_byte = out_byte << 4;
    if (second == -1) out_byte |= 0x02; else out_byte |= (char)second;
    char in_byte = out_byte, tmp_gt;
    tmp_gt = in_byte & 0x03;
    d2 = (tmp_gt == 0x02) ? -1 : (int)tmp_gt;
    in_byte = in_byte >> 4;
    tmp_gt = in_byte & 0x03;
    d1 = (tmp_gt == 0x02) ? -1 : (int)tmp_gt;
}


// ------------------------------------------------------------------------------------------------ per-line summary
struct Summ {
    bool poly_incl = false, poly_all = false, unph_incl = false, unph_all = false;
    long c0 = 0, c1 = 0;
    unsigned Nnm = 0, Nchr = 0;
};

static inline uint32_t popcnt_words(const uint64_t* a, size_t W) {
    uint64_t s = 0;
    for (size_t k = 0; k < W; k++) s += (uint64_t)__builtin_popcountll(a[k]);
    return uint32_t(s);
}

#if defined(__AVX512BW__) && defined(__BMI2__)
#define SWIFTLD_SIMD 1
#include <immintrin.h>
// Regular line fast path: FORMAT is exactly "GT", every sample column is 3 bytes "a?b" with a,b in {0,1,.} and
// ? in {|,/}, columns separated by single tabs (section length 4N-1). Validates every byte; on any deviation
// returns false and the caller uses the generic parser. For these inputs the generic parser's results are:
// ploidy 2, phase = separator, allele = 0/1/-1, enc_dec is the identity, no polyploidy.
// Requires all individuals included (planes are over all samples).
static bool fast_geno_check(const char* p, size_t len, Summ& sm) {
    const __m512i vt = _mm512_set1_epi8('\t'), vb = _mm512_set1_epi8('|'), vs = _mm512_set1_epi8('/');
    const __m512i v0 = _mm512_set1_epi8('0'), v1 = _mm512_set1_epi8('1'), vd = _mm512_set1_epi8('.');
    const uint64_t MT = 0x8888888888888888ull, MS = 0x2222222222222222ull, MA = 0x5555555555555555ull;
    uint64_t c0 = 0, c1 = 0, nsl = 0;
    const size_t nch = (len + 63) / 64;
    for (size_t k = 0; k < nch; k++) {
        const size_t off = k * 64, rem = len - off;
        const uint64_t lm = rem >= 64 ? ~0ull : ((1ull << rem) - 1);
        _mm_prefetch(p + off + 1024, _MM_HINT_T0);   // the line is usually not in cache yet (no newline scan)
        const __m512i v = _mm512_maskz_loadu_epi8(lm, p + off);
        const uint64_t t = _mm512_cmpeq_epi8_mask(v, vt);
        const uint64_t bar = _mm512_cmpeq_epi8_mask(v, vb), sl = _mm512_cmpeq_epi8_mask(v, vs);
        const uint64_t z = _mm512_cmpeq_epi8_mask(v, v0), o = _mm512_cmpeq_epi8_mask(v, v1);
        const uint64_t d = _mm512_cmpeq_epi8_mask(v, vd);
        if (t != (MT & lm)) return false;
        if (((bar | sl) & MS) != (MS & lm)) return false;
        if (((z | o | d) & MA) != (MA & lm)) return false;
        nsl += (uint64_t)__builtin_popcountll(sl & MS);
        c0 += (uint64_t)__builtin_popcountll(z & MA);
        c1 += (uint64_t)__builtin_popcountll(o & MA);
    }
    sm = Summ();
    sm.c0 = long(c0);
    sm.c1 = long(c1);
    sm.Nnm = unsigned(c0 + c1);
    sm.unph_all = sm.unph_incl = (nsl > 0);
    return true;
}
// Planes of a line that fast_geno_check accepted (only for kept sites: most lines fail --maf). pl is zeroed.
static void fast_geno_planes(const char* p, size_t len, bool hap, uint64_t* pl, size_t W) {
    const __m512i v1 = _mm512_set1_epi8('1'), vd = _mm512_set1_epi8('.'), v0 = _mm512_set1_epi8('0');
    const uint64_t MA = 0x5555555555555555ull, M1 = 0x1111111111111111ull, M2 = 0x4444444444444444ull;
    const size_t nplanes = hap ? 2 : 3;
    memset(pl, 0, nplanes * W * sizeof(uint64_t));
    const size_t nch = (len + 63) / 64;
    for (size_t k = 0; k < nch; k++) {
        const size_t off = k * 64, rem = len - off;
        const uint64_t lm = rem >= 64 ? ~0ull : ((1ull << rem) - 1);
        const __m512i v = _mm512_maskz_loadu_epi8(lm, p + off);
        const uint64_t o = _mm512_cmpeq_epi8_mask(v, v1);
        const uint64_t d = _mm512_cmpeq_epi8_mask(v, vd);
        if (hap) {
            const uint64_t A32 = _pext_u64(o, MA), P32 = _pext_u64(~d & lm, MA);
            const unsigned sh = unsigned(32 * (k & 1));
            pl[k >> 1] |= P32 << sh;
            pl[W + (k >> 1)] |= A32 << sh;
        } else {
            const uint64_t z = _mm512_cmpeq_epi8_mask(v, v0);
            const uint64_t nd = ~d & lm;
            const uint64_t V = _pext_u64(nd, M1) & _pext_u64(nd, M2);
            const uint64_t z1 = _pext_u64(z, M1), z2 = _pext_u64(z, M2);
            const unsigned sh = unsigned(16 * (k & 3));
            pl[k >> 2] |= V << sh;
            pl[W + (k >> 2)] |= (V & (z1 | z2)) << sh;
            pl[2 * W + (k >> 2)] |= (V & z1 & z2) << sh;
        }
    }
}
// Generic lines with every individual included and GT the first FORMAT key: one pass over the sample columns
// computes the summary counts and the bit-planes directly (same results as parse_sample into gts, the summary
// loop and the plane loop). Returns false, having changed nothing the caller uses, when the line has fewer than
// N sample columns (lazy retention: the caller's general path handles it).
static bool fused_walk(const char* q, const char* le, bool hap, size_t N, uint64_t* pl, size_t W, Summ& sm) {
    const __m512i vtab = _mm512_set1_epi8('\t');
    auto load_mask = [&](const char* b0) -> uint64_t {
        const size_t rem = size_t(le - b0);
        const __mmask64 lm = rem >= 64 ? ~__mmask64(0) : ((__mmask64(1) << rem) - 1);
        return _mm512_cmpeq_epi8_mask(_mm512_maskz_loadu_epi8(lm, b0), vtab);
    };
    const size_t nplanes = hap ? 2 : 3;
    memset(pl, 0, nplanes * W * sizeof(uint64_t));
    const char* blk = q;
    uint64_t tm = load_mask(blk);
    unsigned pa = 0, ua = 0, nnm = 0, nchr = 0;
    long c0 = 0, c1 = 0;
    uint64_t w0 = 0, w1 = 0, w2 = 0;
    size_t i = 0;
    for (; i < N; i++) {
        while (tm == 0) {
            blk += 64;
            if (blk >= le) break;
            tm = load_mask(blk);
        }
        const char* t = tm ? blk + __builtin_ctzll(tm) : nullptr;
        tm &= tm - 1;
        const char* fend = t ? t : le;
        GT g;
        if (fend - q >= 3 && (q[1] == '/' || q[1] == '|') && q[0] != ':' && q[2] != ':' && (q + 3 == fend || q[3] == ':')) {
            g.ploidy = 2; g.phase = q[1];
            g.a1 = (q[0] != '.') ? (q[0] - '0') : -1;
            g.a2 = (q[2] != '.') ? (q[2] - '0') : -1;
            g.poly = false;
        } else {
            parse_sample(q, fend, 0, g);
        }
        pa |= unsigned(g.poly); ua |= unsigned(g.phase != '|');
        c0 += long(unsigned(g.a1 == 0) + unsigned(g.a2 == 0));
        c1 += long(unsigned(g.a1 == 1) + unsigned(g.a2 == 1));
        nnm += unsigned(g.a1 > -1) + unsigned(g.a2 > -1);
        nchr += g.ploidy;
        int d1, d2;
        enc_dec(g.a1, g.a2, d1, d2);
        if (hap) {
            const unsigned s0 = unsigned((2 * i) & 63);
            w0 |= (uint64_t(d1 >= 0) << s0) | (uint64_t(d2 >= 0) << (s0 + 1));
            w1 |= (uint64_t(d1 > 0) << s0) | (uint64_t(d2 > 0) << (s0 + 1));
            if ((i & 31) == 31) { pl[i >> 5] = w0; pl[W + (i >> 5)] = w1; w0 = w1 = 0; }
        } else {
            const uint64_t v = uint64_t((g.ploidy == 2) & (d1 >= 0) & (d2 >= 0));
            const uint64_t x1 = uint64_t((d1 != d2) | (d1 == 0)), x2 = uint64_t((d1 == d2) & (d1 == 0));
            const unsigned sh = unsigned(i & 63);
            w0 |= v << sh; w1 |= (v & x1) << sh; w2 |= (v & x2) << sh;
            if ((i & 63) == 63) { pl[i >> 6] = w0; pl[W + (i >> 6)] = w1; pl[2 * W + (i >> 6)] = w2; w0 = w1 = w2 = 0; }
        }
        if (!t) { i++; break; }
        q = t + 1;
    }
    if (i < N) return false;
    if (hap) { if (N & 31) { pl[N >> 5] = w0; pl[W + (N >> 5)] = w1; } }
    else if (N & 63) { pl[N >> 6] = w0; pl[W + (N >> 6)] = w1; pl[2 * W + (N >> 6)] = w2; }
    sm = Summ();
    sm.poly_all = sm.poly_incl = pa;
    sm.unph_all = sm.unph_incl = ua;
    sm.c0 = c0; sm.c1 = c1; sm.Nnm = nnm; sm.Nchr = nchr;
    return true;
}
#else
#define SWIFTLD_SIMD 0
#endif

// ------------------------------------------------------------------------------------------------ individuals
static void read_id_file(const string& fn, std::set<string>& out) {
    // filter_individuals_by_keep_list / _exclude_list: while(!eof){getline; ss >> tmp; insert(tmp);}
    FILE* f = fopen(fn.c_str(), "rb");
    if (!f) die("Could not open Individual file:" + fn);
    string content;
    char buf[65536];
    size_t n;
    while ((n = fread(buf, 1, sizeof buf, f)) > 0) content.append(buf, n);
    fclose(f);
    string tmp;
    size_t st = 0;
    for (;;) {
        size_t k = content.find('\n', st);
        sv line = (k == string::npos) ? sv(content).substr(st) : sv(content).substr(st, k - st);
        size_t a = 0;
        auto ws = [](char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r'; };
        while (a < line.size() && ws(line[a])) a++;
        if (a < line.size()) {
            size_t b = a;
            while (b < line.size() && !ws(line[b])) b++;
            tmp = string(line.substr(a, b - a));
        }
        out.insert(tmp);
        if (k == string::npos) break;
        st = k + 1;
    }
}

// ------------------------------------------------------------------------------------------------ output
struct OutBuf {
    char* d = nullptr;
    size_t len = 0, cap = 0;
    bool arena = false;           // d is arena memory (never freed or realloc'ed)
    void reserve(size_t n) {
        if (len + n > cap) {
            size_t nc = std::max(cap * 2, len + n + (1 << 16));
            if (arena) {
                char* nd = (char*)malloc(nc);
                memcpy(nd, d, len);
                d = nd;
                arena = false;
            } else {
                d = (char*)realloc(d, nc);
            }
            cap = nc;
        }
    }
    void free_() { if (!arena) ::free(d); d = nullptr; len = cap = 0; arena = false; }
};

static void write_all(int fd, const char* p, size_t n) {
    while (n) {
        ssize_t w = write(fd, p, n);
        if (w < 0) die("write failed");
        p += w;
        n -= size_t(w);
    }
}

// ------------------------------------------------------------------------------------------------ sites
// Kept sites in fixed blocks, so the parser can append while workers read earlier sites.
struct alignas(64) SiteRec {       // one cache line per site
    int32_t pos = 0;
    uint32_t chrom = 0;           // chromosome id
    uint32_t rank = 0;            // position of this site in its chromosome's index list
    uint32_t c1 = 0, c2 = 0;      // hap: alt count (full sites); geno: sum x, sum x^2 (full sites)
    uint8_t full = 0, pos2len = 0, prelen = 0, prelen_long = 0;
    char pre[24];                 // "CHROM\tPOS\t" when it fits
    char pos2[16];                // "POS\t"
};
// Per-chromosome list of kept-site indices in fixed blocks; the parser appends, workers read entries below the
// published size (release/acquire), so entries never move.
struct ChromList {
    static constexpr size_t B = 4096;
    static constexpr size_t NB = 65536;
    uint32_t** blocks;            // NB block pointers, allocated up front (calloc: untouched until used)
    size_t size = 0;              // written by the parser with __atomic_store_n(release)
    std::atomic<bool> unsorted{false};   // set before pushing a site whose POS is below the previous one's
    ChromList() { blocks = (uint32_t**)calloc(NB, sizeof(uint32_t*)); }
    void push(uint32_t v) {
        if (size % B == 0) {
            if (size / B >= NB) die("too many sites on one chromosome");
            blocks[size / B] = (uint32_t*)malloc(B * sizeof(uint32_t));
        }
        blocks[size / B][size % B] = v;
        __atomic_store_n(&size, size + 1, __ATOMIC_RELEASE);
    }
    uint32_t at(size_t k) const { return blocks[k / B][k % B]; }
    size_t published() const { return __atomic_load_n(&size, __ATOMIC_ACQUIRE); }
};
static_assert(sizeof(SiteRec) == 64, "SiteRec must be one cache line");
struct Sites {
    static constexpr size_t B = 1024;
    struct Block {
        SiteRec rec[B];
        uint64_t* planes;
        const char* lp[B];
        struct RowPart* part[B];   // partially formatted rows (see tile_step), parser-owned until published
        char p2n[B][32];           // "POS\tN\t" with N = all included haplotypes/individuals (complete pairs)
        uint8_t p2nlen[B];
    };
    size_t stride, n = 0, count = 0;
    Block** blocks;
    static constexpr size_t MAXB = size_t(1) << 18;   // 268 M sites
    explicit Sites(size_t st) : stride(st) {
        blocks = (Block**)calloc(MAXB, sizeof(Block*));
    }
    SiteRec& append() {
        if (n % B == 0) {
            if (n / B >= MAXB) die("too many sites");
            void* bm = g_arena.alloc(sizeof(Block));
            Block* b = bm ? new (bm) Block() : new Block();
            b->planes = (uint64_t*)g_arena.alloc((B * stride + 8) * sizeof(uint64_t));
            if (!b->planes) b->planes = (uint64_t*)calloc(B * stride + 8, sizeof(uint64_t));
            blocks[n / B] = b;
        }
        n++;
        return blocks[(n - 1) / B]->rec[(n - 1) % B];
    }
    SiteRec& rec(size_t j) const { return blocks[j / B]->rec[j % B]; }
    const uint64_t* planes(size_t j) const { return blocks[j / B]->planes + (j % B) * stride; }
    uint64_t* planes_mut(size_t j) { return blocks[j / B]->planes + (j % B) * stride; }
    void set_long_prefix(size_t j, const string& s) {
        char* c = (char*)malloc(s.size() + 72);
        uint64_t L = s.size();
        memcpy(c, &L, 8);
        memcpy(c + 8, s.data(), s.size());
        blocks[j / B]->lp[j % B] = c;
    }
    const char* long_prefix(size_t j) const { return blocks[j / B]->lp[j % B] + 8; }
    struct RowPart*& part(size_t j) { return blocks[j / B]->part[j % B]; }
    size_t long_prefix_len(size_t j) const { uint64_t L; memcpy(&L, blocks[j / B]->lp[j % B], 8); return size_t(L); }
};

#if SWIFTFMT_BATCH
// Text of one row's pairs (k2 in [0, m)): "CHROM\tPOS1\t" POS2 "\t" N "\t" r2 [ "\t" D "\t" D' ] "\n".
// A free function with everything passed by value, so the loop keeps its state in registers (inside the
// capturing lambda every captured variable was re-read after each char store). `o` must have room for
// m * (prelen + 100) + 64 bytes. Returns the new end.
template <bool HAP, bool MINR2>
static char* emit_lines(char* o, size_t m, const uint32_t* js, Sites::Block* const* blocks, const char* pre,
                        size_t prelen, bool short_pre, __m256i pre_v, const uint8_t* NT, const uint8_t* NL,
                        uint32_t Hfull, const uint32_t* N_, const double* R2, const double* DD, const double* DP,
                        const uint8_t* GS, const uint16_t* GC, double min_r2) {
    for (size_t k2 = 0; k2 < m; k2++) {
        if (MINR2)
            if ((R2[k2] < min_r2) | (R2[k2] != R2[k2])) continue;
        const uint32_t j = js[k2];
        const Sites::Block* bl = blocks[j / Sites::B];
        const size_t jj = j % Sites::B;
        if (short_pre) _mm256_storeu_si256((__m256i*)o, pre_v); else memcpy(o, pre, prelen);
        o += prelen;
        const uint32_t n = N_[k2];
        if (__builtin_expect(n == Hfull, 1)) {
            // "POS2\tN\t" precomputed per site for complete pairs: one 32-byte copy
            _mm256_storeu_si256((__m256i*)o, _mm256_loadu_si256((const __m256i*)bl->p2n[jj]));
            o += bl->p2nlen[jj];
        } else {
            _mm_storeu_si128((__m128i*)o, _mm_loadu_si128((const __m128i*)bl->rec[jj].pos2));
            o += bl->rec[jj].pos2len;
            // n <= Hfull: "n\t" from the per-run table
            _mm_storeu_si128((__m128i*)o, _mm_loadu_si128((const __m128i*)(NT + 16 * size_t(n))));
            o += NL[n];
        }
        o = swiftfmt::emit16(o, R2[k2], GS + 16 * k2, GC[k2]);
        if (HAP) {
            o = swiftfmt::emit16(o, DD[k2], GS + 16 * (m + k2), GC[m + k2]);
            o = swiftfmt::emit16(o, DP[k2], GS + 16 * (2 * m + k2), GC[2 * m + k2]);
        }
        o[-1] = '\n';
    }
    return o;
}
#endif

// A row whose pairs are formatted ahead of time while the parser waits for decompressed data.
struct RowPart {
    uint32_t k = UINT32_MAX;      // next index in the row's chromosome list to scan (UINT32_MAX: not started)
    bool closed = false;          // the scan has hit its break: the row is complete
    OutBuf buf;                   // text of the pairs formatted so far
};

// Chunks of rows, published by the parser in row order.
struct Chunk {
    uint32_t row_lo = 0, row_hi = 0;
    uint64_t est = 0;             // estimated pairs (chunk sizing and output preallocation only)
    OutBuf* out = nullptr;
    bool done = false;            // guarded by the writer mutex
};
struct ChunkQueue {
    Chunk** slots;
    size_t cap = size_t(1) << 22;
    std::atomic<size_t> published{0}, next_claim{0};
    std::atomic<bool> parse_done{false};
    ChunkQueue() { slots = (Chunk**)calloc(cap, sizeof(Chunk*)); }
    void publish(Chunk* c) {
        size_t k = published.load(std::memory_order_relaxed);
        if (k >= cap) die("too many chunks");
        slots[k] = c;
        published.store(k + 1, std::memory_order_release);
    }
};

static inline uint32_t pop_and(const uint64_t* a, const uint64_t* b, size_t W) {
    uint64_t s = 0;
    for (size_t k = 0; k < W; k++) s += (uint64_t)__builtin_popcountll(a[k] & b[k]);
    return uint32_t(s);
}

// Fused per-pair counts for sites with missing data (one pass over the words instead of 4-9 loops).
// hap: n = |Pi&Pj|, a1 = |Ai&Pj|, a2 = |Aj&Pi|, x22 = |Ai&Aj|
static inline void hap_counts(const uint64_t* Pi, const uint64_t* Ai, const uint64_t* Pj, const uint64_t* Aj, size_t W,
                              uint32_t& n, uint32_t& a1, uint32_t& a2, uint32_t& x22) {
    uint64_t sn = 0, s1 = 0, s2 = 0, sx = 0;
    for (size_t k = 0; k < W; k++) {
        sn += (uint64_t)__builtin_popcountll(Pi[k] & Pj[k]);
        s1 += (uint64_t)__builtin_popcountll(Ai[k] & Pj[k]);
        s2 += (uint64_t)__builtin_popcountll(Aj[k] & Pi[k]);
        sx += (uint64_t)__builtin_popcountll(Ai[k] & Aj[k]);
    }
    n = uint32_t(sn); a1 = uint32_t(s1); a2 = uint32_t(s2); x22 = uint32_t(sx);
}
// geno, dosage x = a + b with b subset of a subset of V. Sxy = |ai&aj| + |ai&bj| + |bi&aj| + |bi&bj|
// = |ai&aj| + |(ai&bj)|(bi&aj)| + 2|bi&bj|, because (ai&bj)&(bi&aj) = bi&bj.
static inline uint32_t geno_xy(const uint64_t* ai, const uint64_t* bi, const uint64_t* aj, const uint64_t* bj, size_t W) {
    uint64_t s = 0;
    for (size_t k = 0; k < W; k++) {
        s += (uint64_t)__builtin_popcountll(ai[k] & aj[k]);
        s += (uint64_t)__builtin_popcountll((ai[k] & bj[k]) | (bi[k] & aj[k]));
        s += 2 * (uint64_t)__builtin_popcountll(bi[k] & bj[k]);
    }
    return uint32_t(s);
}
static inline void geno_counts(const uint64_t* Vi, const uint64_t* ai, const uint64_t* bi, const uint64_t* Vj,
                               const uint64_t* aj, const uint64_t* bj, size_t W, uint32_t& n, uint32_t& sX,
                               uint32_t& sX2, uint32_t& sY, uint32_t& sY2, uint32_t& sXY) {
    uint64_t sn = 0, aiv = 0, biv = 0, ajv = 0, bjv = 0, xy = 0;
    for (size_t k = 0; k < W; k++) {
        sn += (uint64_t)__builtin_popcountll(Vi[k] & Vj[k]);
        aiv += (uint64_t)__builtin_popcountll(ai[k] & Vj[k]);
        biv += (uint64_t)__builtin_popcountll(bi[k] & Vj[k]);
        ajv += (uint64_t)__builtin_popcountll(aj[k] & Vi[k]);
        bjv += (uint64_t)__builtin_popcountll(bj[k] & Vi[k]);
        xy += (uint64_t)__builtin_popcountll(ai[k] & aj[k]);
        xy += (uint64_t)__builtin_popcountll((ai[k] & bj[k]) | (bi[k] & aj[k]));
        xy += 2 * (uint64_t)__builtin_popcountll(bi[k] & bj[k]);
    }
    n = uint32_t(sn);
    sX = uint32_t(aiv + biv); sX2 = uint32_t(aiv + 3 * biv);
    sY = uint32_t(ajv + bjv); sY2 = uint32_t(ajv + 3 * bjv);
    sXY = uint32_t(xy);
}

#if defined(__AVX512F__) && defined(__AVX512VPOPCNTDQ__) && defined(__AVX512VL__)
#define SWIFTLD_VCOUNT 1
// Lane sums of eight vectors: result lane t = sum of the 8 qwords of v[t] (exact integer adds).
static inline __m512i hsum8x8(const __m512i v[8]) {
    const __m512i a0 = _mm512_add_epi64(_mm512_unpacklo_epi64(v[0], v[1]), _mm512_unpackhi_epi64(v[0], v[1]));
    const __m512i a1 = _mm512_add_epi64(_mm512_unpacklo_epi64(v[2], v[3]), _mm512_unpackhi_epi64(v[2], v[3]));
    const __m512i a2 = _mm512_add_epi64(_mm512_unpacklo_epi64(v[4], v[5]), _mm512_unpackhi_epi64(v[4], v[5]));
    const __m512i a3 = _mm512_add_epi64(_mm512_unpacklo_epi64(v[6], v[7]), _mm512_unpackhi_epi64(v[6], v[7]));
    const __m512i b0 = _mm512_add_epi64(_mm512_shuffle_i64x2(a0, a1, 0x88), _mm512_shuffle_i64x2(a0, a1, 0xDD));
    const __m512i b1 = _mm512_add_epi64(_mm512_shuffle_i64x2(a2, a3, 0x88), _mm512_shuffle_i64x2(a2, a3, 0xDD));
    return _mm512_add_epi64(_mm512_shuffle_i64x2(b0, b1, 0x88), _mm512_shuffle_i64x2(b0, b1, 0xDD));
}
// Counts for up to 8 pairs at once (row site i against sites with plane pointers pj[t]), all sites complete.
// hap: out[t] = |Ai & Aj|.  geno: out[t] = |ai&aj| + |(ai&bj)|(bi&aj)| + 2|bi&bj|  (planes: V, a, b).
template <bool HAP>
static inline void vcount8(const uint64_t* pi, const uint64_t* const* pj, size_t nb, size_t W, uint32_t* out) {
    // pj[t] for t >= nb point at a zero plane block (contributes 0); the t loops have constant trip counts so
    // the accumulators stay in registers
    __m512i acc0 = _mm512_setzero_si512(), acc1 = acc0, acc2 = acc0, acc3 = acc0, acc4 = acc0, acc5 = acc0,
            acc6 = acc0, acc7 = acc0;
    const uint64_t *q0 = pj[0], *q1 = pj[1], *q2 = pj[2], *q3 = pj[3], *q4 = pj[4], *q5 = pj[5], *q6 = pj[6], *q7 = pj[7];
    for (size_t w = 0; w < W; w += 8) {
        const __mmask8 lm = (W - w >= 8) ? 0xff : __mmask8((1u << (W - w)) - 1);
        if (HAP) {
            const __m512i ai = _mm512_maskz_loadu_epi64(lm, pi + W + w);
#define VC_H(acc, q) acc = _mm512_add_epi64(acc, _mm512_popcnt_epi64(_mm512_and_si512(ai, _mm512_maskz_loadu_epi64(lm, q + W + w))))
            VC_H(acc0, q0); VC_H(acc1, q1); VC_H(acc2, q2); VC_H(acc3, q3);
            VC_H(acc4, q4); VC_H(acc5, q5); VC_H(acc6, q6); VC_H(acc7, q7);
#undef VC_H
        } else {
            const __m512i ai = _mm512_maskz_loadu_epi64(lm, pi + W + w), bi = _mm512_maskz_loadu_epi64(lm, pi + 2 * W + w);
#define VC_G(acc, q) { const __m512i aj = _mm512_maskz_loadu_epi64(lm, q + W + w), bj = _mm512_maskz_loadu_epi64(lm, q + 2 * W + w); \
                const __m512i p1 = _mm512_popcnt_epi64(_mm512_and_si512(ai, aj)); \
                const __m512i p2 = _mm512_popcnt_epi64(_mm512_or_si512(_mm512_and_si512(ai, bj), _mm512_and_si512(bi, aj))); \
                const __m512i p3 = _mm512_popcnt_epi64(_mm512_and_si512(bi, bj)); \
                acc = _mm512_add_epi64(acc, _mm512_add_epi64(_mm512_add_epi64(p1, p2), _mm512_add_epi64(p3, p3))); }
            VC_G(acc0, q0); VC_G(acc1, q1); VC_G(acc2, q2); VC_G(acc3, q3);
            VC_G(acc4, q4); VC_G(acc5, q5); VC_G(acc6, q6); VC_G(acc7, q7);
#undef VC_G
        }
    }
    const __m512i v[8] = {acc0, acc1, acc2, acc3, acc4, acc5, acc6, acc7};
    const __m512i r = hsum8x8(v);
    _mm256_mask_storeu_epi32(out, __mmask8((1u << nb) - 1), _mm512_cvtepi64_epi32(r));
}
// Counts for 8 pairs that may have missing data, packed four 16-bit fields per 64-bit lane so that one transposed
// horizontal sum serves four counts (exact while every count < 65536: hap H <= 65535, geno 4 NI <= 65535).
// hap: n = |Pi&Pj|, a1 = |Ai&Pj|, a2 = |Aj&Pi|, x22 = |Ai&Aj|  -> N_, A_, B_, X_
// geno: n = |Vi&Vj|, sX = |ai&Vj| + |bi&Vj|, sX2 = |ai&Vj| + 3|bi&Vj|, sY, sY2 likewise, sXY  -> N_, A_, B_, X_, Y_, Z_
template <bool HAP>
static inline void vcount8_part(const uint64_t* pi, const uint64_t* const* pj, size_t nb, size_t W, size_t k2,
                                uint32_t* N_, uint32_t* A_, uint32_t* B_, uint32_t* X_, uint32_t* Y_, uint32_t* Z_) {
    __m512i p0[8], p1[8];
    for (int t = 0; t < 8; t++) { p0[t] = _mm512_setzero_si512(); p1[t] = _mm512_setzero_si512(); }
    for (size_t w = 0; w < W; w += 8) {
        const __mmask8 lm = (W - w >= 8) ? 0xff : __mmask8((1u << (W - w)) - 1);
        if (HAP) {
            const __m512i Pi = _mm512_maskz_loadu_epi64(lm, pi + w), Ai = _mm512_maskz_loadu_epi64(lm, pi + W + w);
#pragma GCC unroll 8
            for (int t = 0; t < 8; t++) {
                const __m512i Pj = _mm512_maskz_loadu_epi64(lm, pj[t] + w), Aj = _mm512_maskz_loadu_epi64(lm, pj[t] + W + w);
                __m512i v = _mm512_popcnt_epi64(_mm512_and_si512(Pi, Pj));
                v = _mm512_or_si512(v, _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(Ai, Pj)), 16));
                v = _mm512_or_si512(v, _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(Aj, Pi)), 32));
                v = _mm512_or_si512(v, _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(Ai, Aj)), 48));
                p0[t] = _mm512_add_epi64(p0[t], v);
            }
        } else {
            const __m512i Vi = _mm512_maskz_loadu_epi64(lm, pi + w), ai = _mm512_maskz_loadu_epi64(lm, pi + W + w),
                          bi = _mm512_maskz_loadu_epi64(lm, pi + 2 * W + w);
#pragma GCC unroll 8
            for (int t = 0; t < 8; t++) {
                const __m512i Vj = _mm512_maskz_loadu_epi64(lm, pj[t] + w), aj = _mm512_maskz_loadu_epi64(lm, pj[t] + W + w),
                              bj = _mm512_maskz_loadu_epi64(lm, pj[t] + 2 * W + w);
                __m512i v = _mm512_popcnt_epi64(_mm512_and_si512(Vi, Vj));
                v = _mm512_or_si512(v, _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(ai, Vj)), 16));
                v = _mm512_or_si512(v, _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(bi, Vj)), 32));
                v = _mm512_or_si512(v, _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(aj, Vi)), 48));
                p0[t] = _mm512_add_epi64(p0[t], v);
                const __m512i xy = _mm512_add_epi64(
                    _mm512_add_epi64(_mm512_popcnt_epi64(_mm512_and_si512(ai, aj)),
                                     _mm512_popcnt_epi64(_mm512_or_si512(_mm512_and_si512(ai, bj), _mm512_and_si512(bi, aj)))),
                    _mm512_slli_epi64(_mm512_popcnt_epi64(_mm512_and_si512(bi, bj)), 1));
                __m512i u = _mm512_popcnt_epi64(_mm512_and_si512(bj, Vi));
                u = _mm512_or_si512(u, _mm512_slli_epi64(xy, 16));
                p1[t] = _mm512_add_epi64(p1[t], u);
            }
        }
    }
    const __mmask8 om = __mmask8((1u << nb) - 1);
    const __m512i f = _mm512_set1_epi64(0xffff);
    const __m512i r0 = hsum8x8(p0);
    const __m512i c0 = _mm512_and_si512(r0, f), c1 = _mm512_and_si512(_mm512_srli_epi64(r0, 16), f),
                  c2 = _mm512_and_si512(_mm512_srli_epi64(r0, 32), f), c3 = _mm512_srli_epi64(r0, 48);
    _mm256_mask_storeu_epi32(N_ + k2, om, _mm512_cvtepi64_epi32(c0));
    if (HAP) {
        _mm256_mask_storeu_epi32(A_ + k2, om, _mm512_cvtepi64_epi32(c1));
        _mm256_mask_storeu_epi32(B_ + k2, om, _mm512_cvtepi64_epi32(c2));
        _mm256_mask_storeu_epi32(X_ + k2, om, _mm512_cvtepi64_epi32(c3));
    } else {
        const __m512i r1 = hsum8x8(p1);
        const __m512i bjv = _mm512_and_si512(r1, f), sxy = _mm512_srli_epi64(r1, 16);
        // c1 = aiv, c2 = biv, c3 = ajv
        _mm256_mask_storeu_epi32(A_ + k2, om, _mm512_cvtepi64_epi32(_mm512_add_epi64(c1, c2)));
        _mm256_mask_storeu_epi32(B_ + k2, om, _mm512_cvtepi64_epi32(_mm512_add_epi64(c1, _mm512_add_epi64(c2, _mm512_add_epi64(c2, c2)))));
        _mm256_mask_storeu_epi32(X_ + k2, om, _mm512_cvtepi64_epi32(_mm512_add_epi64(c3, bjv)));
        _mm256_mask_storeu_epi32(Y_ + k2, om, _mm512_cvtepi64_epi32(_mm512_add_epi64(c3, _mm512_add_epi64(bjv, _mm512_add_epi64(bjv, bjv)))));
        _mm256_mask_storeu_epi32(Z_ + k2, om, _mm512_cvtepi64_epi32(sxy));
    }
}
#else
#define SWIFTLD_VCOUNT 0
#endif

#if SWIFTLD_VCOUNT
// Rows with missing data, kept out of the hot row function (its code layout matters for the complete-data rows).
template <bool HAP>
__attribute__((noinline)) static size_t partial_counts(const uint64_t* pi, const uint32_t* js, size_t m, const Sites& SS,
                                                       size_t W, const uint64_t* ZP, uint32_t Hfull, uint32_t* N_,
                                                       uint32_t* A_, uint32_t* B_, uint32_t* X_, uint32_t* Y_, uint32_t* Z_) {
    const uint64_t* pj8[8];
    for (size_t k2 = 0; k2 < m; k2 += 8) {
        const size_t nb = std::min<size_t>(8, m - k2);
        for (size_t t = 0; t < 8; t++) pj8[t] = t < nb ? SS.planes(js[k2 + t]) : ZP;
        vcount8_part<HAP>(pi, pj8, nb, W, k2, N_, A_, B_, X_, Y_, Z_);
    }
    size_t np = 0;
    for (size_t k2 = 0; k2 < m; k2++) np += (N_[k2] != Hfull);
    return np;
}
#endif
#if defined(__AVX512F__) && defined(__AVX512VL__)
// calc_hap_r2 for general rows, 8 pairs at a time: the five quotients by cc = n come from t2 when every lane has
// H - n < T2D (t2[d * T2S + k] = k / (H - d), the same IEEE quotients), else from divisions; then the
// reference's expressions in order.
__attribute__((noinline)) static void arith_general_hap(size_t m, const uint32_t* N_, const uint32_t* A_, const uint32_t* B_,
                                                        const uint32_t* X_, uint32_t Hfull, const double* t2, size_t T2D,
                                                        size_t T2S, double* R2, double* DD, double* DP) {
    // general rows (some n < H), 8 pairs at a time: the five quotients by cc = n come from T2 when
    // every lane has H - n < T2D, else from divisions; then calc_hap_r2's expressions in order
    
    const __m512d zero = _mm512_setzero_pd();
    const __m256i vH = _mm256_set1_epi32(int(Hfull)), vD = _mm256_set1_epi32(int(T2D)),
                  vS = _mm256_set1_epi32(int(T2S));
    for (size_t k2 = 0; k2 < m; k2 += 8) {
        const __mmask8 lm = (m - k2 >= 8) ? 0xff : __mmask8((1u << (m - k2)) - 1);
        const __m256i n = _mm256_maskz_loadu_epi32(lm, N_ + k2), A = _mm256_maskz_loadu_epi32(lm, A_ + k2);
        const __m256i B = _mm256_maskz_loadu_epi32(lm, B_ + k2), X = _mm256_maskz_loadu_epi32(lm, X_ + k2);
        const __m256i i21 = _mm256_sub_epi32(A, X), i12 = _mm256_sub_epi32(B, X);
        const __m256i i11 = _mm256_sub_epi32(_mm256_sub_epi32(_mm256_sub_epi32(n, X), i21), i12);
        const __m256i k1 = _mm256_add_epi32(i11, i12), k2v = _mm256_add_epi32(i21, X);
        const __m256i k3 = _mm256_add_epi32(i11, i21), k4 = _mm256_add_epi32(i12, X);
        const __m256i d = _mm256_sub_epi32(vH, n);
        __m512d rel, p1, p2, q1, q2;
        if ((_mm256_mask_cmplt_epu32_mask(lm, d, vD)) == lm) {
            const __m256i base = _mm256_mullo_epi32(d, vS);
            rel = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, i11), t2, 8);
            p1 = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, k1), t2, 8);
            p2 = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, k2v), t2, 8);
            q1 = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, k3), t2, 8);
            q2 = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, k4), t2, 8);
        } else {
            // (x11 + x12) etc. are exact integer sums, so converting the int sums is identical
            const __m512d cc = _mm512_cvtepi32_pd(n);
            rel = _mm512_div_pd(_mm512_cvtepi32_pd(i11), cc);
            p1 = _mm512_div_pd(_mm512_cvtepi32_pd(k1), cc);
            p2 = _mm512_div_pd(_mm512_cvtepi32_pd(k2v), cc);
            q1 = _mm512_div_pd(_mm512_cvtepi32_pd(k3), cc);
            q2 = _mm512_div_pd(_mm512_cvtepi32_pd(k4), cc);
        }
        const __m512d p1q1 = _mm512_mul_pd(p1, q1);
        const __m512d D = _mm512_sub_pd(rel, p1q1);
        const __m512d m1 = _mm512_min_pd(_mm512_mul_pd(p2, q2), p1q1);
        const __m512d m2 = _mm512_min_pd(_mm512_mul_pd(p2, q1), _mm512_mul_pd(p1, q2));
        const __m512d Dmax = _mm512_mask_blend_pd(_mm512_cmp_pd_mask(D, zero, _CMP_LT_OQ), m2, m1);
        _mm512_mask_storeu_pd(DD + k2, lm, D);
        _mm512_mask_storeu_pd(DP + k2, lm, _mm512_div_pd(D, Dmax));
        const __m512d var1 = _mm512_sub_pd(p1, _mm512_mul_pd(p1, p1));
        const __m512d var2 = _mm512_sub_pd(q1, _mm512_mul_pd(q1, q1));
        const __m512d cov = _mm512_sub_pd(rel, p1q1);
        _mm512_mask_storeu_pd(R2 + k2, lm, _mm512_div_pd(_mm512_mul_pd(cov, cov), _mm512_mul_pd(var1, var2)));
    }
}
// calc_geno_r2 for general rows (same scheme)
__attribute__((noinline)) static void arith_general_geno(size_t m, const uint32_t* N_, const uint32_t* A_, const uint32_t* B_,
                                                         const uint32_t* X_, const uint32_t* Y_, const uint32_t* Z_,
                                                         uint32_t Hfull, const double* t2, size_t T2D, size_t T2S, double* R2) {
    
    const __m512d zero = _mm512_setzero_pd();
    const __m256i vH = _mm256_set1_epi32(int(Hfull)), vD = _mm256_set1_epi32(int(T2D)),
                  vS = _mm256_set1_epi32(int(T2S));
    for (size_t k2 = 0; k2 < m; k2 += 8) {
        const __mmask8 lm = (m - k2 >= 8) ? 0xff : __mmask8((1u << (m - k2)) - 1);
        const __m256i n = _mm256_maskz_loadu_epi32(lm, N_ + k2);
        const __m256i sx = _mm256_maskz_loadu_epi32(lm, A_ + k2), sx2 = _mm256_maskz_loadu_epi32(lm, B_ + k2);
        const __m256i sy = _mm256_maskz_loadu_epi32(lm, X_ + k2), sy2 = _mm256_maskz_loadu_epi32(lm, Y_ + k2);
        const __m256i sxy = _mm256_maskz_loadu_epi32(lm, Z_ + k2);
        const __m256i d = _mm256_sub_epi32(vH, n);
        __m512d X, X2, Y, Y2, XY;
        if ((_mm256_mask_cmplt_epu32_mask(lm, d, vD)) == lm) {
            const __m256i base = _mm256_mullo_epi32(d, vS);
            X = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, sx), t2, 8);
            X2 = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, sx2), t2, 8);
            Y = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, sy), t2, 8);
            Y2 = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, sy2), t2, 8);
            XY = _mm512_mask_i32gather_pd(zero, lm, _mm256_add_epi32(base, sxy), t2, 8);
        } else {
            const __m512d ic = _mm512_cvtepi32_pd(n);
            X = _mm512_div_pd(_mm512_cvtepu32_pd(sx), ic);
            X2 = _mm512_div_pd(_mm512_cvtepu32_pd(sx2), ic);
            Y = _mm512_div_pd(_mm512_cvtepu32_pd(sy), ic);
            Y2 = _mm512_div_pd(_mm512_cvtepu32_pd(sy2), ic);
            XY = _mm512_div_pd(_mm512_cvtepu32_pd(sxy), ic);
        }
        const __m512d var1 = _mm512_sub_pd(X2, _mm512_mul_pd(X, X));
        const __m512d var2 = _mm512_sub_pd(Y2, _mm512_mul_pd(Y, Y));
        const __m512d cov = _mm512_sub_pd(XY, _mm512_mul_pd(X, Y));
        _mm512_mask_storeu_pd(R2 + k2, lm, _mm512_div_pd(_mm512_mul_pd(cov, cov), _mm512_mul_pd(var1, var2)));
    }
}
#endif

static int n_threads_env() {
    const char* e = getenv("SWIFTLD_THREADS");
    int t = e ? atoi(e) : 0;
    if (t <= 0) {
        cpu_set_t cs;
        CPU_ZERO(&cs);
        if (sched_getaffinity(0, sizeof cs, &cs) == 0) t = CPU_COUNT(&cs);
        if (t <= 0) t = (int)sysconf(_SC_NPROCESSORS_ONLN);
        if (t <= 0) t = 1;
    }
    return std::min(t, 256);
}

#include <chrono>
static double now_s() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv) {
    const bool prof = getenv("SWIFTLD_PROFILE") != nullptr;
    double t_start = now_s();
    g_argv = argv;
    Params P = parse_args(argc, argv);
    const bool allow_stream = n_threads_env() >= 2 && getenv("SWIFTLD_NOSTREAM") == nullptr;
    // One malloc arena for all threads. Otherwise each thread's first malloc creates an arena (a 128 MB
    // PROT_NONE reservation, then two munmap trims), which needs the mmap lock for writing; on gzip input the
    // inflater thread did this at its start and waited behind the main thread's page population (several ms).
    mallopt(M_ARENA_MAX, 1);
    g_arena.init(size_t(1) << 30);   // address space only; pages are faulted on use
    Input in = read_input(P, allow_stream);
    double t_read = now_s();
    LineReader lr{in.data, in.data + in.len, in.st, nullptr};
    // With page recycling the inflater's later output lands in recycled pages; pages populated ahead of it
    // would be replaced (populated, then freed by mremap), so do not prefault at all, also not while the
    // header is being read.
    if (in.st && swiftld_block_hook) lr.no_prefault = true;

    // ---- header (vcf_file::read_header + header::parse_header)
    vector<string> indv;
    bool has_genotypes = false;
    {
        const char *s, *e;
        while (lr.next(s, e)) {
            size_t n = size_t(e - s);
            if (n > 2 && s[0] == '#') {
                if (s[1] == '#') continue;  // meta line
                vector<sv> f;
                split_sv(sv(s, n), '\t', f);
                has_genotypes = false;
                for (size_t c = 0; c < f.size(); c++) {
                    if (c == 8) { if (f[c] == "FORMAT") has_genotypes = true; }
                    else if (c > 8) indv.push_back(string(f[c]));
                }
                break;
            }
            break;  // a non-header line ends the header and is consumed
        }
    }
    const size_t N = indv.size();
    vector<char> include(N, 1);
    // filter_individuals(..., keep_then_exclude=true)
    if (!P.keep_files.empty()) {
        std::set<string> keep;
        for (auto& f : P.keep_files) read_id_file(f, keep);
        for (size_t i = 0; i < N; i++) if (include[i] && !keep.count(indv[i])) include[i] = 0;
    }
    if (!P.remove_files.empty()) {
        std::set<string> rem;
        for (auto& f : P.remove_files) read_id_file(f, rem);
        for (size_t i = 0; i < N; i++) if (include[i] && rem.count(indv[i])) include[i] = 0;
    }
    vector<int> incl_idx;
    for (size_t i = 0; i < N; i++) if (include[i]) incl_idx.push_back(int(i));
    const size_t NI = incl_idx.size();
    if (!has_genotypes || NI == 0) die("Require Genotypes in VCF file in order to output LD Statistics.");

    // ---- output target (created before the site loop, as the reference does)
    int out_fd = 1;
    if (!P.stream_out) {
        string fn = P.out + (P.hap ? ".hap.ld" : ".geno.ld");
        out_fd = open(fn.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (out_fd < 0) die("Could not open LD Output File: " + fn, 3);
    }

    // ---- site store: stable addresses, because pair workers read sites while the parser appends them
    const bool hap = P.hap;
    const size_t H = hap ? 2 * NI : NI;  // bit slots per site
    const size_t W = (H + 63) / 64;
    const size_t nplanes = hap ? 2 : 3;
    Sites SS(nplanes * W);
    std::unordered_map<string, int> chrom_ids;
    vector<string> chrom_names;
    vector<int> s_chrom, s_pos;                 // parser-owned copies (row completion test, chunk sizing)
    vector<ChromList*> clists;                  // shared with workers; reserved so it never reallocates
    clists.reserve(1 << 20);
    vector<vector<uint32_t>> by_chrom;          // parser-owned, for chunk-size estimates only
    vector<uint32_t> rank_in_chrom;
    vector<char> chrom_sorted;

    // retained-field state (std::getline into member strings keeps old values when the stream is exhausted)
    sv R_ID, R_REF, R_ALT, R_FILTER;
    sv R_FORMAT;                    // last fully parsed FORMAT_str
    // last fully parsed genotype columns: either a whole line (lazy) or a materialised vector
    const char* R_geno_line = nullptr; const char* R_geno_end = nullptr;
    bool R_geno_materialised = true;
    vector<sv> R_geno(N, sv());
    string thin_chrom; int thin_pos = 0;
    long N_kept_entries = 0;
    const bool freq_filter = !((P.min_maf <= 0.0) && (P.max_maf >= 1.0) && (P.min_site_call_rate <= 0));
    vector<GT> gts(N);
    const bool all_included = (NI == N);
    size_t n_fast_lines = 0, n_generic_lines = 0;
    vector<uint64_t> scratch(nplanes * W + 8, 0);
    vector<sv> tmpv;
    vector<sv> cur_cols;

    // ---- pair scheduling. Rows (ui) are cut into chunks of ~3000 pairs as soon as their pair list is final:
    // the reference's scan of row ui stops at the first SNP-window or bp-window break, so once a parsed site
    // proves that break, later sites cannot add pairs to the row. The parser enumerates each finished row
    // (exact j list) and publishes chunks in row order; workers format chunks and write them in order.
    const unsigned skip = (unsigned)std::max(1, P.ld_snp_window_min);
    const int snp_window = P.ld_snp_window_size, bp_min = P.ld_bp_window_min, bp_max = P.ld_bp_window_size;
    const double min_r2 = P.min_r2;
    ChunkQueue Q;
    const int T_threads_hint = n_threads_env();
    size_t next_row = 0;
    Chunk* cur = nullptr;
    const size_t TARGET = 8000;
    auto row_complete = [&](size_t r, size_t L) -> bool {
        if (int64_t(L) - int64_t(r) >= int64_t(snp_window)) return true;          // every later j breaks
        if (s_chrom[L] != s_chrom[r]) return false;
        if (uint64_t(L) < uint64_t(r) + skip) return false;
        if (int64_t(L) - int64_t(r) > int64_t(snp_window)) return false;
        const int d = int(unsigned(s_pos[L]) - unsigned(s_pos[r]));
        return d >= bp_min && d > bp_max;                                          // site L is a bp break
    };
    uint64_t est_acc = 0;
    auto advance_rows = [&](bool final) {
        const size_t count = s_pos.size();
        if (count < 2) return;
        const size_t L = count - 1;
        while (next_row < count - 1 && (final || row_complete(next_row, L))) {
            if (!cur) { cur = new Chunk(); cur->row_lo = uint32_t(next_row); }
            // chunk size estimate only (never affects output): same-chromosome candidates inside the SNP
            // window, cut at the bp window by binary search while the chromosome is sorted so far
            {
                const size_t r = next_row, cid = size_t(s_chrom[r]);
                const vector<uint32_t>& cl = by_chrom[cid];
                uint64_t est = 0;
                if (snp_window >= 1) {
                    size_t a = size_t(std::lower_bound(cl.begin() + rank_in_chrom[r], cl.end(), uint64_t(r) + skip,
                                           [](uint32_t x, uint64_t y) { return uint64_t(x) < y; }) - cl.begin());
                    size_t b = size_t(std::upper_bound(cl.begin() + a, cl.end(), uint64_t(r) + uint64_t(snp_window),
                                           [](uint64_t y, uint32_t x) { return y < uint64_t(x); }) - cl.begin());
                    if (chrom_sorted[cid] && b > a) {
                        const int64_t lim = int64_t(s_pos[r]) + int64_t(bp_max);
                        b = size_t(std::upper_bound(cl.begin() + a, cl.begin() + b, lim,
                                       [&](int64_t y, uint32_t x) { return y < int64_t(s_pos[x]); }) - cl.begin());
                    }
                    est = b > a ? b - a : 0;
                }
                est_acc += est + 8;
                cur->est += est;
            }
            next_row++;
            cur->row_hi = uint32_t(next_row);
            if (est_acc >= TARGET) { Q.publish(cur); cur = nullptr; est_acc = 0; }
        }
    };

    // ---- output machinery
    const bool gate_writes = P.stream_out && in.st != nullptr;   // stdout cannot be rewritten after a re-exec
    std::mutex mu;
    size_t next_write = 0;
    bool writer_active = false, header_written = false;
    size_t out_bytes = 0, falloc_end = 0;   // bytes written to the output file; preallocated up to falloc_end
    uint64_t est_written = 0;               // estimated pairs of the chunks written so far
    vector<OutBuf*> pool;
    char Hstr[32] = {0};
    size_t Hlen = size_t(swiftfmt::fmt_u32(Hstr, uint32_t(hap ? H : NI)) - Hstr);
    Hstr[Hlen++] = '\t';
    const uint32_t Hfull = uint32_t(hap ? H : NI);
    // "n\t" for every possible pair count n <= Hfull (16-byte slots) and its length
    vector<uint8_t> Ntext(16 * (size_t(Hfull) + 1), 0), Nlen(size_t(Hfull) + 1, 0);
    for (uint32_t n = 0; n <= Hfull; n++) {
        char* q = swiftfmt::fmt_u32((char*)&Ntext[16 * n], n);
        *q++ = '\t';
        Nlen[n] = uint8_t(q - (char*)&Ntext[16 * n]);
    }
    const uint8_t* const NT = Ntext.data();
    const uint8_t* const NL = Nlen.data();
    // Quotients k / Hfull for every count k a pair of complete sites can produce (hap: k <= H; geno: k <= 4 NI).
    // For such pairs the reference divides integer-valued doubles by double(Hfull); the table holds the same
    // IEEE quotients, so lookups are bitwise identical to the divisions they replace.
    vector<uint64_t> zero_planes(nplanes * W + 8, 0);   // padding target for vcount8
    const uint64_t* const ZP = zero_planes.data();
    vector<double> Tq(size_t(hap ? H : 4 * NI) + 1);
    { const double den = double(int(Hfull)); for (size_t k = 0; k < Tq.size(); k++) Tq[k] = double(int(k)) / den; }
    const double* const TQ = Tq.data();
    // Pairs with missing calls divide by n = Hfull - d. For d < T2D the quotients k / n come from a second table
    // (rows of length T2S, at most 1 MB), built on first use: T2[d * T2S + k] = double(k) / double(Hfull - d).
    const size_t T2S = size_t(hap ? Hfull : 4 * Hfull) + 1;
    const size_t T2D = std::min<size_t>({size_t(64), size_t(Hfull) + 1, (size_t(1) << 17) / T2S});
    vector<double> T2v;
    std::once_flag t2_once;
    auto T2 = [&]() -> const double* {
        std::call_once(t2_once, [&]() {
            T2v.resize(T2D * T2S);
            for (size_t d = 0; d < T2D; d++) {
                const double den = double(int(Hfull - d));
                for (size_t k = 0; k < T2S; k++) T2v[d * T2S + k] = double(int(k)) / den;
            }
        });
        return T2v.data();
    };
    const char* hdr = hap ? "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime\n" : "CHR\tPOS1\tPOS2\tN_INDV\tR^2\n";
    std::atomic<long> ns_fmt{0}, ns_write{0};
    auto get_buf = [&]() -> OutBuf* {
        {
            std::lock_guard<std::mutex> g(mu);
            if (!pool.empty()) { OutBuf* b = pool.back(); pool.pop_back(); b->len = 0; return b; }
        }
        OutBuf* b = new OutBuf();
        if (char* m = (char*)g_arena.alloc(size_t(1) << 20)) {
            b->d = m; b->cap = size_t(1) << 20; b->arena = true;
            for (size_t o = 0; o < b->cap; o += 4096) ((volatile char*)m)[o] = 0;   // pre-fault once; recycled
        } else {
            b->reserve(1 << 20);
            memset(b->d, 0, b->cap);  // pre-fault once; the buffer is recycled
        }
        return b;
    };
    // writes every finished chunk that is next in order; one writer at a time; caller holds no lock
    auto flush_ready = [&]() {
        std::unique_lock<std::mutex> lk(mu);
        if (writer_active) return;
        if (gate_writes && !g_inflate_ok.load(std::memory_order_acquire)) return;
        writer_active = true;
        for (;;) {
            const size_t pubd = Q.published.load(std::memory_order_acquire);
            if (next_write >= pubd) break;
            Chunk* ch = Q.slots[next_write];
            if (!ch->done) break;
            const size_t hlen = header_written ? 0 : strlen(hdr);
            // Reserve file blocks ahead (beyond EOF, file size unchanged): buffered writes into allocated blocks
            // skip ext4's per-page delayed-allocation work (~8% faster in a 150 MB test). 32 MB steps while the
            // input is being parsed; once every chunk is published, only bytes already formatted (done chunks in
            // order from this one) are reserved, so the reservation normally ends inside the file and the final
            // ftruncate that would release an unused tail (0.25-0.7 ms, serial at the end) is not needed.
            size_t want = 0;
            if (!P.stream_out && falloc_end != SIZE_MAX && out_bytes + hlen + ch->out->len > falloc_end) {
                want = size_t(32) << 20;
                if (Q.parse_done.load(std::memory_order_acquire)) {
                    // every chunk is published: bytes of formatted chunks plus 0.9 x (estimated pairs of the
                    // others) x (bytes per pair so far). An overestimate (e.g. --min-r2) only means the final
                    // ftruncate is still needed.
                    uint64_t done_b = 0, done_e = 0, rem_e = 0;
                    for (size_t c = next_write; c < pubd; c++) {
                        const Chunk* x = Q.slots[c];
                        if (x->done) { done_b += x->out->len; done_e += x->est; } else rem_e += x->est;
                    }
                    const double bpp = double(out_bytes + done_b) / double(std::max<uint64_t>(1, est_written + done_e));
                    const uint64_t total = out_bytes + hlen + done_b + uint64_t(0.9 * bpp * double(rem_e));
                    want = total > falloc_end ? std::min<uint64_t>(want, total - falloc_end) : 0;
                    if (want < (size_t(2) << 20)) want = 0;   // not worth a call; these bytes are written unreserved
                }
            }
            lk.unlock();
            double tw0 = prof ? now_s() : 0;
            if (!header_written) { write_all(out_fd, hdr, hlen); header_written = true; out_bytes += hlen; }
            if (want) {
                if (fallocate(out_fd, FALLOC_FL_KEEP_SIZE, off_t(falloc_end), off_t(want)) == 0) falloc_end += want;
                else falloc_end = SIZE_MAX;   // not supported here: stop trying
            }
            write_all(out_fd, ch->out->d, ch->out->len);
            out_bytes += ch->out->len;
            est_written += ch->est;
            if (prof) ns_write += long((now_s() - tw0) * 1e9);
            lk.lock();
            pool.push_back(ch->out);
            ch->out = nullptr;
            next_write++;
        }
        writer_active = false;
    };

    // per-thread scratch for the row pipeline
    struct RowScratch {
        vector<uint32_t> cn, ca, cb, cx, cy, cz;
        vector<double> vr2, vD, vDp;
        vector<uint8_t> g_s16;   // 16-byte field sources: [0,m) r2, [m,2m) D, [2m,3m) D'
        vector<uint16_t> g_c16;
        vector<uint32_t> jv;
    };
    // The reference's scan of row ui over same-chromosome sites (others only `continue`), from chromosome-list
    // index k (UINT32_MAX: start of the row) up to `ncl` entries. Returns the next index; *closed is set when the
    // scan hits its break. Sites beyond the parser's frontier are never dereferenced: a scan stops at ncl, and a
    // complete row breaks on the index test (SNP window) or at a parsed bp-window break before reaching them.
    auto scan_row = [&](size_t ui, uint32_t kstart, size_t ncl, vector<uint32_t>& jv, bool* closed) -> uint32_t {
        const SiteRec& r0 = SS.rec(ui);
        const int posi = r0.pos;
        const ChromList& cl = *clists[r0.chrom];
        size_t k = kstart;
        if (k == UINT32_MAX) {
            // first entry with index >= ui + skip (entries are increasing site indices)
            k = r0.rank + 1;
            const uint64_t lo = uint64_t(ui) + skip;
            if (k < ncl && uint64_t(cl.at(k)) < lo) {
                size_t a = k, b = ncl;
                while (a < b) { size_t mid = (a + b) / 2; if (uint64_t(cl.at(mid)) < lo) a = mid + 1; else b = mid; }
                k = a;
            }
        }
        *closed = false;
        // Fast path: the remaining entries are consecutive site indices and the chromosome's positions so far are
        // non-decreasing. Then the SNP-window break is at a computable index, d = POS_j - POS_i is monotone, the
        // `continue`d entries (d < bp_min) form a prefix and the bp break is the first d > bp_max: binary searches
        // give exactly the entries and the break index of the loop below.
        if (k < ncl && !cl.unsorted.load(std::memory_order_acquire) &&
            size_t(cl.at(ncl - 1)) - size_t(cl.at(k)) == ncl - 1 - k && SS.rec(cl.at(k)).pos >= posi) {
            const size_t base = cl.at(k);   // site index of entry k; entry kk is site base + (kk - k)
            const int64_t lim = int64_t(ui) + int64_t(snp_window) - int64_t(base);   // entries kk - k <= lim pass
            const size_t ksnp = lim < 0 ? k : size_t(std::min<int64_t>(int64_t(ncl - k), lim + 1)) + k;
            auto dpos = [&](size_t kk) -> int64_t { return int64_t(SS.rec(base + (kk - k)).pos) - int64_t(posi); };
            // (the reference computes d as int(unsigned difference); identical while all d fit in an int)
            if (ksnp > k && dpos(ksnp - 1) > int64_t(INT_MAX)) goto general;
            auto first_ge = [&](int64_t t) -> size_t {   // first kk in [k, ksnp) with d >= t (d is non-decreasing)
                size_t a = k, b = ksnp;
                while (a < b) { const size_t mid = (a + b) / 2; if (dpos(mid) < t) a = mid + 1; else b = mid; }
                return a;
            };
            const size_t kb = first_ge(int64_t(bp_min));
            const size_t kc = first_ge(std::max<int64_t>(int64_t(bp_min), int64_t(bp_max) + 1));
            for (size_t kk = kb; kk < kc; kk++) jv.push_back(uint32_t(base + (kk - k)));
            if (kc < ksnp) { *closed = true; return uint32_t(kc); }
            if (ksnp < ncl) { *closed = true; return uint32_t(ksnp); }
            return uint32_t(ncl);
        }
    general:
        for (; k < ncl; k++) {
            const size_t uj = cl.at(k);
            if (int(uj - ui) > snp_window) { *closed = true; break; }
            const int d = int(unsigned(SS.rec(uj).pos) - unsigned(posi));
            if (d < bp_min) continue;
            if (d > bp_max) { *closed = true; break; }
            jv.push_back(uint32_t(uj));
        }
        return uint32_t(k);
    };
    // counts -> reference arithmetic -> text for the pairs (ui, js[0..m)), appended to ob
    auto format_row = [&](size_t ui, const uint32_t* js, size_t m, OutBuf& ob, RowScratch& R) {
        vector<uint32_t>& cn = R.cn; vector<uint32_t>& ca = R.ca; vector<uint32_t>& cb = R.cb;
        vector<uint32_t>& cx = R.cx; vector<uint32_t>& cy = R.cy; vector<uint32_t>& cz = R.cz;
        vector<double>& vr2 = R.vr2; vector<double>& vD = R.vD; vector<double>& vDp = R.vDp;
        vector<uint8_t>& g_s16 = R.g_s16; vector<uint16_t>& g_c16 = R.g_c16;
            if (cn.size() < m) {
                size_t nm = std::max(m, cn.size() * 2);
                cn.resize(nm); ca.resize(nm); cb.resize(nm); cx.resize(nm); cy.resize(nm); cz.resize(nm);
                vr2.resize(nm); vD.resize(nm); vDp.resize(nm);
                g_s16.resize(16 * (3 * nm + 8)); g_c16.resize(3 * nm + 8);
            }
            const SiteRec& ri = SS.rec(ui);
            const uint64_t* pi = SS.planes(ui);
            const bool fi = ri.full;
            uint32_t* __restrict N_ = cn.data();
            uint32_t* __restrict A_ = ca.data();
            uint32_t* __restrict B_ = cb.data();
            uint32_t* __restrict X_ = cx.data();
            uint32_t* __restrict Y_ = cy.data();
            uint32_t* __restrict Z_ = cz.data();
            double* __restrict R2 = vr2.data();
            double* __restrict DD = vD.data();
            double* __restrict DP = vDp.data();
            if (hap) {
                // counts: N_=n, A_=a1 (alt at i among pairwise present), B_=a2, X_=x22
                const uint32_t ci = ri.c1;
                size_t n_partial = 0;
                bool vec_done = false;
#if SWIFTLD_VCOUNT
                if (fi) {
                    // pairs of complete sites: x22 = |Ai & Aj| for 8 pairs at a time (one pass: site record,
                    // plane pointer and count); a row with any incomplete site falls back to the loop below
                    const uint64_t* pj8[8];
                    bool allf = true;
                    for (size_t k2 = 0; k2 < m && allf; k2 += 8) {
                        const size_t nb = std::min<size_t>(8, m - k2);
                        for (size_t t = 0; t < 8; t++) {
                            if (t < nb) {
                                const uint32_t j = js[k2 + t];
                                const Sites::Block* bl = SS.blocks[j / Sites::B];
                                const SiteRec& rj = bl->rec[j % Sites::B];
                                allf &= bool(rj.full);
                                B_[k2 + t] = rj.c1;
                                pj8[t] = bl->planes + (j % Sites::B) * SS.stride;
                            } else pj8[t] = ZP;
                        }
                        if (allf) vcount8<true>(pi, pj8, nb, W, X_ + k2);
                    }
                    if (allf) {
                        for (size_t k2 = 0; k2 < m; k2++) { N_[k2] = Hfull; A_[k2] = ci; }
                        vec_done = true;
                    }
                }
#endif
#if SWIFTLD_VCOUNT
                if (!vec_done && Hfull <= 65535) {
                    // some pair has missing data: all four counts for 8 pairs at a time (packed 16-bit fields)
                    n_partial = partial_counts<true>(pi, js, m, SS, W, ZP, Hfull, N_, A_, B_, X_, nullptr, nullptr);
                    vec_done = true;
                }
#endif
                if (!vec_done)
                for (size_t k2 = 0; k2 < m; k2++) {
                    const size_t uj = js[k2];
                    const SiteRec& rj = SS.rec(uj);
                    const uint64_t* pj = SS.planes(uj);
                    if (fi && rj.full) { N_[k2] = Hfull; A_[k2] = ci; B_[k2] = rj.c1; X_[k2] = pop_and(pi + W, pj + W, W); }
                    else { hap_counts(pi, pi + W, pj, pj + W, W, N_[k2], A_[k2], B_[k2], X_[k2]); n_partial++; }
                }
#if defined(__AVX512F__) && defined(__AVX512VL__)
                if (n_partial == 0) {
                    // every pair has n = H: x11 / cc, (x11+x12) / cc, ... are table lookups (same IEEE quotients);
                    // p1 = (H - ci) / H and p2 = ci / H are row constants. Then calc_hap_r2's expressions, in order.
                    const double p1s = TQ[Hfull - ci], p2s = TQ[ci];
                    const __m512d p1 = _mm512_set1_pd(p1s), p2 = _mm512_set1_pd(p2s);
                    const __m512d var1 = _mm512_set1_pd(p1s - p1s * p1s), zero = _mm512_setzero_pd();
                    const __m256i vhc = _mm256_set1_epi32(int(Hfull - ci)), vH = _mm256_set1_epi32(int(Hfull));
                    for (size_t k2 = 0; k2 < m; k2 += 8) {
                        const __mmask8 lm = (m - k2 >= 8) ? 0xff : __mmask8((1u << (m - k2)) - 1);
                        const __m256i X = _mm256_maskz_loadu_epi32(lm, X_ + k2), B = _mm256_maskz_loadu_epi32(lm, B_ + k2);
                        const __m256i i11 = _mm256_add_epi32(_mm256_sub_epi32(vhc, B), X);   // H - ci - cj + x22
                        const __m512d rel = _mm512_mask_i32gather_pd(zero, lm, i11, TQ, 8);
                        const __m512d q1 = _mm512_mask_i32gather_pd(zero, lm, _mm256_sub_epi32(vH, B), TQ, 8);
                        const __m512d q2 = _mm512_mask_i32gather_pd(zero, lm, B, TQ, 8);
                        const __m512d p1q1 = _mm512_mul_pd(p1, q1);
                        const __m512d D = _mm512_sub_pd(rel, p1q1);
                        // std::min(a, b) == (b < a) ? b : a == MINPD(b, a)
                        const __m512d m1 = _mm512_min_pd(_mm512_mul_pd(p2, q2), p1q1);
                        const __m512d m2 = _mm512_min_pd(_mm512_mul_pd(p2, q1), _mm512_mul_pd(p1, q2));
                        const __m512d Dmax = _mm512_mask_blend_pd(_mm512_cmp_pd_mask(D, zero, _CMP_LT_OQ), m2, m1);
                        _mm512_mask_storeu_pd(DD + k2, lm, D);
                        _mm512_mask_storeu_pd(DP + k2, lm, _mm512_div_pd(D, Dmax));
                        const __m512d var2 = _mm512_sub_pd(q1, _mm512_mul_pd(q1, q1));
                        const __m512d cov = _mm512_sub_pd(rel, p1q1);
                        _mm512_mask_storeu_pd(R2 + k2, lm, _mm512_div_pd(_mm512_mul_pd(cov, cov), _mm512_mul_pd(var1, var2)));
                    }
                } else
#endif
#if defined(__AVX512F__) && defined(__AVX512VL__)
                if (n_partial != 0) {
                    arith_general_hap(m, N_, A_, B_, X_, Hfull, T2(), T2D, T2S, R2, DD, DP);
                } else
#endif
                // calc_hap_r2, same double expression order
                for (size_t k2 = 0; k2 < m; k2++) {
                    const int ix22 = int(X_[k2]), ix21 = int(A_[k2]) - ix22, ix12 = int(B_[k2]) - ix22;
                    const int ix11 = int(N_[k2]) - ix22 - ix21 - ix12;
                    const double x11 = ix11, x12 = ix12, x21 = ix21, x22d = ix22;
                    const double cc = double(int(N_[k2]));
                    const double rel_x11 = x11 / cc;
                    const double p1 = (x11 + x12) / cc;
                    const double p2 = (x21 + x22d) / cc;
                    const double q1 = (x11 + x21) / cc;
                    const double q2 = (x12 + x22d) / cc;
                    const double D = rel_x11 - p1 * q1;
                    const double m1 = std::min(p1 * q1, p2 * q2), m2 = std::min(p1 * q2, p2 * q1);
                    const double Dmax = (D < 0) ? m1 : m2;
                    DD[k2] = D;
                    DP[k2] = D / Dmax;
                    // X = (x11+x12)/cc == p1, X2 == X, Y == q1, Y2 == Y, XY == rel_x11 (same operands, same rounding)
                    const double var1 = p1 - p1 * p1, var2 = q1 - q1 * q1, cov12 = rel_x11 - p1 * q1;
                    R2[k2] = cov12 * cov12 / (var1 * var2);
                }
            } else {
                // counts: N_=n, A_=sX, B_=sX2, X_=sY, Y_=sY2, Z_=sXY
                const uint64_t *Vi = pi, *ai = pi + W, *bi = pi + 2 * W;
                size_t n_partial = 0;
                bool vec_done = false;
#if SWIFTLD_VCOUNT
                if (fi) {
                    const uint64_t* pj8[8];
                    bool allf = true;
                    for (size_t k2 = 0; k2 < m && allf; k2 += 8) {
                        const size_t nb = std::min<size_t>(8, m - k2);
                        for (size_t t = 0; t < 8; t++) {
                            if (t < nb) {
                                const uint32_t j = js[k2 + t];
                                const Sites::Block* bl = SS.blocks[j / Sites::B];
                                const SiteRec& rj = bl->rec[j % Sites::B];
                                allf &= bool(rj.full);
                                X_[k2 + t] = rj.c1; Y_[k2 + t] = rj.c2;
                                pj8[t] = bl->planes + (j % Sites::B) * SS.stride;
                            } else pj8[t] = ZP;
                        }
                        if (allf) vcount8<false>(pi, pj8, nb, W, Z_ + k2);
                    }
                    if (allf) {
                        for (size_t k2 = 0; k2 < m; k2++) { N_[k2] = Hfull; A_[k2] = ri.c1; B_[k2] = ri.c2; }
                        vec_done = true;
                    }
                }
#endif
#if SWIFTLD_VCOUNT
                if (!vec_done && 4 * size_t(Hfull) <= 65535) {
                    n_partial = partial_counts<false>(pi, js, m, SS, W, ZP, Hfull, N_, A_, B_, X_, Y_, Z_);
                    vec_done = true;
                }
#endif
                if (!vec_done)
                for (size_t k2 = 0; k2 < m; k2++) {
                    const size_t uj = js[k2];
                    const SiteRec& rj = SS.rec(uj);
                    const uint64_t* pj = SS.planes(uj);
                    const uint64_t *Vj = pj, *aj = pj + W, *bj = pj + 2 * W;
                    if (fi && rj.full) {
                        Z_[k2] = geno_xy(ai, bi, aj, bj, W);
                        N_[k2] = Hfull; A_[k2] = ri.c1; B_[k2] = ri.c2; X_[k2] = rj.c1; Y_[k2] = rj.c2;
                    } else {
                        geno_counts(Vi, ai, bi, Vj, aj, bj, W, N_[k2], A_[k2], B_[k2], X_[k2], Y_[k2], Z_[k2]);
                        n_partial++;
                    }
                }
#if defined(__AVX512F__) && defined(__AVX512VL__)
                if (n_partial == 0) {
                    // every pair has n = NI: the five quotients are table lookups (same IEEE quotients);
                    // X, X2 and var1 are row constants. Then calc_geno_r2's expressions, in order.
                    const double Xs = TQ[ri.c1], X2s = TQ[ri.c2];
                    const __m512d X = _mm512_set1_pd(Xs), var1 = _mm512_set1_pd(X2s - Xs * Xs), zero = _mm512_setzero_pd();
                    for (size_t k2 = 0; k2 < m; k2 += 8) {
                        const __mmask8 lm = (m - k2 >= 8) ? 0xff : __mmask8((1u << (m - k2)) - 1);
                        const __m512d Y = _mm512_mask_i32gather_pd(zero, lm, _mm256_maskz_loadu_epi32(lm, X_ + k2), TQ, 8);
                        const __m512d Y2 = _mm512_mask_i32gather_pd(zero, lm, _mm256_maskz_loadu_epi32(lm, Y_ + k2), TQ, 8);
                        const __m512d XY = _mm512_mask_i32gather_pd(zero, lm, _mm256_maskz_loadu_epi32(lm, Z_ + k2), TQ, 8);
                        const __m512d var2 = _mm512_sub_pd(Y2, _mm512_mul_pd(Y, Y));
                        const __m512d cov = _mm512_sub_pd(XY, _mm512_mul_pd(X, Y));
                        _mm512_mask_storeu_pd(R2 + k2, lm, _mm512_div_pd(_mm512_mul_pd(cov, cov), _mm512_mul_pd(var1, var2)));
                    }
                } else
#endif
#if defined(__AVX512F__) && defined(__AVX512VL__)
                if (n_partial != 0) {
                    arith_general_geno(m, N_, A_, B_, X_, Y_, Z_, Hfull, T2(), T2D, T2S, R2);
                } else
#endif
                // calc_geno_r2
                for (size_t k2 = 0; k2 < m; k2++) {
                    const double ic = double(int(N_[k2]));
                    const double X = double(A_[k2]) / ic, X2 = double(B_[k2]) / ic;
                    const double Y = double(X_[k2]) / ic, Y2 = double(Y_[k2]) / ic, XY = double(Z_[k2]) / ic;
                    const double var1 = X2 - X * X, var2 = Y2 - Y * Y, cov12 = XY - X * Y;
                    R2[k2] = cov12 * cov12 / (var1 * var2);
                }
            }
#if SWIFTFMT_BATCH
            // digits and layout codes of all values of this row, 8 lanes at a time
            uint8_t* GS = g_s16.data();
            uint16_t* GC = g_c16.data();
            swiftfmt::g6_batch16(R2, m, GS, GC);
            if (hap) {
                swiftfmt::g6_batch16(DD, m, GS + 16 * m, GC + m);
                swiftfmt::g6_batch16(DP, m, GS + 32 * m, GC + 2 * m);
            }
            // output: row prefix and the N field are row constants (kept in registers); every field is
            // stored with one 16- or 32-byte store and the pointer advanced by its length
            const bool short_pre = !ri.prelen_long;
            const char* pre = short_pre ? ri.pre : SS.long_prefix(ui);
            const size_t prelen = short_pre ? ri.prelen : SS.long_prefix_len(ui);
            const __m256i pre_v = _mm256_loadu_si256((const __m256i*)ri.pre);   // 24-byte slot + 8 bytes of the record
            ob.reserve(m * (std::max<size_t>(prelen, 32) + 100) + 64);
            char* o = ob.d + ob.len;
            if (hap) o = (min_r2 > 0) ? emit_lines<true, true>(o, m, js, SS.blocks, pre, prelen, short_pre, pre_v, NT, NL, Hfull, N_, R2, DD, DP, GS, GC, min_r2)
                                      : emit_lines<true, false>(o, m, js, SS.blocks, pre, prelen, short_pre, pre_v, NT, NL, Hfull, N_, R2, DD, DP, GS, GC, min_r2);
            else o = (min_r2 > 0) ? emit_lines<false, true>(o, m, js, SS.blocks, pre, prelen, short_pre, pre_v, NT, NL, Hfull, N_, R2, DD, DP, GS, GC, min_r2)
                                  : emit_lines<false, false>(o, m, js, SS.blocks, pre, prelen, short_pre, pre_v, NT, NL, Hfull, N_, R2, DD, DP, GS, GC, min_r2);
            ob.len = size_t(o - ob.d);
#else
            const char* pre = ri.prelen_long ? SS.long_prefix(ui) : ri.pre;
            const size_t prelen = ri.prelen_long ? SS.long_prefix_len(ui) : ri.prelen;
            for (size_t k2 = 0; k2 < m; k2++) {
                const double r2 = R2[k2];
                if (min_r2 > 0)
                    if ((r2 < min_r2) | (r2 != r2)) continue;
                const SiteRec& rj = SS.rec(js[k2]);
                ob.reserve(prelen + 128);
                char* o = ob.d + ob.len;
                memcpy(o, pre, prelen);
                o += prelen;
                memcpy(o, rj.pos2, 16);
                o += rj.pos2len;
                const uint32_t n = N_[k2];
                if (n == Hfull) { memcpy(o, Hstr, 16); o += Hlen; }
                else { o = swiftfmt::fmt_u32(o, n); *o++ = '\t'; }
                o = swiftfmt::fmt_g(o, r2);
                if (hap) {
                    *o++ = '\t';
                    o = swiftfmt::fmt_g(o, DD[k2]);
                    *o++ = '\t';
                    o = swiftfmt::fmt_g(o, DP[k2]);
                }
                *o++ = '\n';
                ob.len = size_t(o - ob.d);
            }
#endif
    };
    auto process_chunk = [&](size_t c) {
        double tc0 = prof ? now_s() : 0;
        Chunk* ch = Q.slots[c];
        OutBuf* obp = get_buf();
        OutBuf& ob = *obp;
        thread_local RowScratch R;
        for (size_t ui = ch->row_lo; ui < ch->row_hi; ui++) {
            RowPart* part = SS.part(ui);
            uint32_t k0 = UINT32_MAX;
            if (part) {
                if (part->buf.len) {
                    ob.reserve(part->buf.len + 64);
                    memcpy(ob.d + ob.len, part->buf.d, part->buf.len);
                    ob.len += part->buf.len;
                }
                const bool cl_ = part->closed;
                k0 = part->k;
                part->buf.free_();
                delete part;
                SS.part(ui) = nullptr;
                if (cl_) continue;
            }
            R.jv.clear();
            bool closed;
            scan_row(ui, k0, clists[SS.rec(ui).chrom]->published(), R.jv, &closed);
            if (!R.jv.empty()) format_row(ui, R.jv.data(), R.jv.size(), ob, R);
        }
        if (prof) ns_fmt += long((now_s() - tc0) * 1e9);
        {
            std::lock_guard<std::mutex> g(mu);
            ch->out = obp;
            ch->done = true;
        }
        flush_ready();
    };
    // claim the next published chunk, if any
    auto try_one = [&]() -> bool {
        size_t c = Q.next_claim.load(std::memory_order_relaxed);
        if (c >= Q.published.load(std::memory_order_acquire)) return false;
        if (!Q.next_claim.compare_exchange_strong(c, c + 1)) return true;  // lost the race; caller retries
        process_chunk(c);
        return true;
    };
    // ---- genotype pre-parsing (plain VCF, all individuals included). The sample section of a generic line
    // (FORMAT starting with GT, not the regular 4N-1 layout) costs microseconds per line; while the parser works
    // on line L, an otherwise idle worker runs fused_walk on lines L+1 .. L+RING into a ring of slots. The
    // parser uses a slot only if it holds exactly its current line (sequence number and line start); otherwise
    // it computes the line itself. fused_walk is a pure function of the line bytes, so the result is the same.
    struct PreSlot {
        std::atomic<uint64_t> seq{0};   // line ordinal + 1 when the slot holds that line's result
        const char* ls = nullptr;
        bool ok = false;
        Summ sm;
    };
    const size_t PRE_RING = std::max<size_t>(16, std::min<size_t>(256, (size_t(1) << 17) / (nplanes * W + 8)));
    const bool pre_enabled = SWIFTLD_SIMD && !in.st && all_included && T_threads_hint >= 2;   // uses fused_walk
    vector<PreSlot> pre_slots(pre_enabled ? PRE_RING : 0);
    vector<uint64_t> pre_planes(pre_enabled ? PRE_RING * (nplanes * W + 8) : 0);
    std::atomic<uint64_t> parser_line{0};    // lines the parser has finished
    std::atomic<bool> pre_busy{false};
    LineReader pre_lr{nullptr, nullptr, nullptr, nullptr};
    uint64_t pre_line = 0;                   // next line ordinal for the pre-parser (owner: pre_busy holder)
    std::atomic<bool> pre_ready{false};
    vector<sv> pre_tmp;
    auto preparse_step = [&]() -> bool {
        if (!pre_ready.load(std::memory_order_acquire)) return false;
        bool expect = false;
        if (!pre_busy.compare_exchange_strong(expect, true)) return false;
        bool did = false;
        for (int rep = 0; rep < 4; rep++) {
            if (pre_line >= parser_line.load(std::memory_order_acquire) + PRE_RING - 1) break;   // ring full
            const char *ls, *le;
            if (!pre_lr.next(ls, le)) { pre_line = UINT64_MAX / 2; break; }
            const uint64_t L = pre_line++;
            if (L < parser_line.load(std::memory_order_acquire)) continue;   // the parser is already past it
            // same tokenization as the parse loop: the first 9 fields, then the sample section [p, le)
            const char* fs8 = nullptr; const char* fe8 = nullptr;
            const char* p = ls; int k = 0; bool geno_cols = false;
            while (k < 9) {
                const char* t = (const char*)memchr(p, '\t', size_t(le - p));
                if (!t) { if (k == 8) { fs8 = p; fe8 = le; } k++; p = le; break; }
                if (k == 8) { fs8 = p; fe8 = t; }
                k++; p = t + 1;
                if (k == 9) geno_cols = true;
            }
            if (!geno_cols) continue;
            const sv FORMAT(fs8, size_t(fe8 - fs8));
            int gt = -1;
            if (FORMAT == "GT") gt = 0;
            else if (!FORMAT.empty()) {
                split_sv(FORMAT, ':', pre_tmp);
                for (size_t c = 0; c < pre_tmp.size(); c++) if (pre_tmp[c] == "GT") gt = int(c);
            }
            if (gt != 0) continue;
            if (FORMAT.size() == 2 && size_t(le - p) == 4 * N - 1) continue;   // the parser's AVX-512 path takes it
#if SWIFTLD_SIMD
            PreSlot& sl = pre_slots[L % PRE_RING];
            uint64_t* pl = pre_planes.data() + (L % PRE_RING) * (nplanes * W + 8);
            sl.ls = ls;
            sl.ok = fused_walk(p, le, hap, N, pl, W, sl.sm);
            sl.seq.store(L + 1, std::memory_order_release);
            did = true;
#endif
        }
        pre_busy.store(false, std::memory_order_release);
        return did;
    };
    auto worker_loop = [&]() {
        int idle = 0;
        for (;;) {
            if (try_one()) { idle = 0; continue; }
            if (pre_enabled && !Q.parse_done.load(std::memory_order_acquire) && preparse_step()) { idle = 0; continue; }
            if (Q.parse_done.load(std::memory_order_acquire) &&
                Q.next_claim.load() >= Q.published.load(std::memory_order_acquire)) break;
            if (++idle < 2000) __builtin_ia32_pause();
            else { sched_yield(); if (idle > 20000) { struct timespec ts = {0, 50000}; nanosleep(&ts, nullptr); } }
        }
    };

    // ---- threads: [inflater -> worker] (streaming gz), main = parser (helps while waiting) -> worker, extras
    const int T = n_threads_env();
    vector<std::thread> extra;
    if (in.st) {
        g_after_inflate_fn = [&]() { flush_ready(); worker_loop(); };
        g_after_inflate_ready.store(true, std::memory_order_release);
    }
    {
        int n_extra = T - 1 - (in.st ? 1 : 0);
        for (int t = 0; t < n_extra; t++) extra.emplace_back(worker_loop);
    }
    // Rows that are not complete yet can still be formatted up to the parsed frontier while the parser waits for
    // decompressed data (all-pairs files have no complete row until the end of the input). The text goes into a
    // per-row buffer; the worker that later processes the row copies it and continues the scan where it stopped.
    size_t tile_row = 0;
    RowScratch tileR;
    // extension step in chromosome-list entries; SWIFTLD_TILE_STEP is a test hook (1 = extend at every chance)
    const size_t tile_step_entries = getenv("SWIFTLD_TILE_STEP") ? size_t(std::max(1, atoi(getenv("SWIFTLD_TILE_STEP")))) : 512;
    auto tile_step = [&]() -> bool {
        const size_t count = s_pos.size();
        if (count < 2 || next_row + 1 >= count) return false;
        size_t budget = 4096;
        bool did = false;
        size_t r = std::max(tile_row, next_row);
        for (size_t visited = 0; visited < count && budget > 0; visited++) {
            if (r + 1 >= count) { r = next_row; if (did) break; }
            RowPart*& pp = SS.part(r);
            if (!pp) pp = new RowPart();
            if (!pp->closed) {
                const size_t ncl = clists[SS.rec(r).chrom]->size;
                const size_t from = pp->k == UINT32_MAX ? SS.rec(r).rank + 1 : pp->k;
                if (ncl >= from + tile_step_entries) {   // large steps only (per-call overhead, copy later)
                    tileR.jv.clear();
                    bool closed;
                    pp->k = scan_row(r, pp->k, ncl, tileR.jv, &closed);
                    pp->closed = closed;
                    if (!tileR.jv.empty()) {
                        format_row(r, tileR.jv.data(), tileR.jv.size(), pp->buf, tileR);
                        budget = tileR.jv.size() >= budget ? 0 : budget - tileR.jv.size();
                        did = true;
                    }
                }
            }
            r++;
        }
        tile_row = r;
        return did;
    };
    lr.idle = [&]() -> bool { return try_one() || tile_step(); };

    // page recycling (see Stream): publish the lowest input offset this loop may still read
    Stream* const rst = (in.st && swiftld_block_hook) ? in.st : nullptr;
    const char* R_geno_min = nullptr;   // lowest address in R_geno when it is materialised
    auto lowest_retained = [&](const char* lo) -> const char* {
        auto f = [&](sv v) { if (!v.empty() && v.data() < lo) lo = v.data(); };
        f(R_ID); f(R_REF); f(R_ALT); f(R_FILTER); f(R_FORMAT);
        if (R_geno_materialised) { if (R_geno_min && R_geno_min < lo) lo = R_geno_min; }
        else if (R_geno_line && R_geno_line < lo) lo = R_geno_line;
        return lo;
    };
    if (pre_enabled) {
        pre_lr.p = lr.p;
        pre_lr.end = lr.end;
        pre_ready.store(true, std::memory_order_release);
    }
    uint64_t line_no = 0;
    const char *ls, *le;
    // regular lines (FORMAT GT, every sample included, 3-byte calls) skip the newline search (LineReader::next_pred)
    // (plain VCF only: with gzip input the parser reads lines the inflater has just written, which are in cache,
    // so the search is cheap there, and the predicted path measured slower on r_all_hap, exp 52)
    const size_t n4 = (all_included && N > 0 && !in.st) ? 4 * N - 1 : 0;
    bool prov = false;   // the current line's end is predicted and not yet verified
    while (lr.next_pred(ls, le, n4, prov)) {
        const uint64_t cur_line = line_no++;
        if (pre_enabled) parser_line.store(cur_line, std::memory_order_release);
        if (rst) rst->safe.store(size_t(lowest_retained(ls) - rst->buf), std::memory_order_release);
        if (ls == le) continue;  // entry::apply_filters: empty line never passes, not counted
        // tokenize first 9 fields
        const char* fs[9]; const char* fe[9]; int k = 0;
        const char* p = ls;
        bool geno_cols = false;
        auto tokenize = [&]() {
            k = 0; p = ls; geno_cols = false;
            while (k < 9) {
                const char* t = (const char*)memchr(p, '\t', size_t(le - p));
                if (!t) { fs[k] = p; fe[k] = le; k++; p = le; break; }
                fs[k] = p; fe[k] = t; k++; p = t + 1;
                if (k == 9) geno_cols = true;
            }
        };
        tokenize();
        // a predicted line end is accepted only for a regular line that passes the regular-line check (which
        // admits only genotype characters, so no newline); otherwise the true end is found and, if it moved,
        // the fixed fields are tokenized again
        bool pre_ok = false;
        Summ sm_pre;
        if (prov) {
            prov = false;
#if SWIFTLD_SIMD
            if (geno_cols && size_t(fe[8] - fs[8]) == 2 && fs[8][0] == 'G' && fs[8][1] == 'T' && size_t(le - p) == n4 &&
                fast_geno_check(p, n4, sm_pre))
                pre_ok = true;
            else
#endif
            {
                const char* old = le;
                lr.settle(ls, le);
                if (le != old) tokenize();
            }
        }
        sv CHROM(fs[0], size_t(fe[0] - fs[0]));
        sv ID = (k >= 2) ? sv(fs[1], size_t(fe[1] - fs[1])) : R_ID;
        int POS;
        { char b[32]; size_t m = std::min<size_t>(ID.size(), 31); memcpy(b, ID.data(), m); b[m] = 0; POS = atoi(b); }
        if (k >= 3) ID = sv(fs[2], size_t(fe[2] - fs[2]));
        R_ID = ID;
        sv REF = (k >= 4) ? sv(fs[3], size_t(fe[3] - fs[3])) : R_REF; R_REF = REF;
        sv ALT = (k >= 5) ? sv(fs[4], size_t(fe[4] - fs[4])) : R_ALT; R_ALT = ALT;
        sv FILTER = (k >= 7) ? sv(fs[6], size_t(fe[6] - fs[6])) : R_FILTER; R_FILTER = FILTER;

        // ALT alleles (vcf_entry::set_ALT + entry::add_ALT_allele: "." tokens are dropped)
        split_sv(ALT, ',', tmpv);
        size_t n_alt = 0;
        bool alt_len_differs = false;
        for (auto& t : tmpv) if (t != ".") { n_alt++; if (t.size() != REF.size()) alt_len_differs = true; }
        const size_t N_alleles = n_alt + 1;

        // 1. --remove-indels (filter_sites_by_allele_type)
        if (P.remove_indels) {
            bool is_indel = (REF.size() != 1) || alt_len_differs;
            if (is_indel) continue;
        }
        // 2. --remove-filtered-all (filter_sites_by_filter_status)
        if (P.remove_filtered_all) {
            if (FILTER != ".") {
                split_sv(FILTER, ';', tmpv);
                sv mn = tmpv[0];
                for (auto& t : tmpv) if (t < mn) mn = t;
                if (mn != "PASS") continue;
            }
        }
        // 3. --chr (filter_sites_by_chromosome)
        if (!P.chrs.empty() && !P.chrs.count(string(CHROM))) continue;
        // 4. min_alleles = max_alleles = 2 (set by --hap-r2/--geno-r2)
        if (N_alleles != 2) continue;

        // ---- genotype stage: parsed once per line (vcftools re-parses up to three times, identically)
        int GT_idx = -2;
        bool fast = false, fused = false;
        Summ sm;
        auto ensure_parsed = [&]() {
            if (GT_idx != -2) return;
            sv FORMAT = (k >= 9) ? sv(fs[8], size_t(fe[8] - fs[8])) : R_FORMAT;
            R_FORMAT = FORMAT;
            GT_idx = -1;
            if (FORMAT == "GT") GT_idx = 0;
            else if (!FORMAT.empty()) {
                split_sv(FORMAT, ':', tmpv);
                for (size_t c = 0; c < tmpv.size(); c++) if (tmpv[c] == "GT") GT_idx = int(c);
            }
#if SWIFTLD_SIMD
            if (geno_cols && all_included && GT_idx == 0 && FORMAT.size() == 2 && size_t(le - p) == 4 * N - 1) {
                if (pre_ok ? (sm = sm_pre, true) : fast_geno_check(p, size_t(le - p), sm)) {
                    sm.Nchr = unsigned(2 * N);
                    fast = true;
                    n_fast_lines++;
                    R_geno_materialised = false;
                    R_geno_line = p;
                    R_geno_end = le;
                    return;
                }
            }
#endif
#if SWIFTLD_SIMD
            if (pre_enabled && geno_cols && GT_idx == 0) {
                PreSlot& sl = pre_slots[cur_line % PRE_RING];
                if (sl.seq.load(std::memory_order_acquire) == cur_line + 1 && sl.ls == ls) {
                    if (sl.ok) {
                        sm = sl.sm;
                        memcpy(scratch.data(), pre_planes.data() + (cur_line % PRE_RING) * (nplanes * W + 8),
                               nplanes * W * sizeof(uint64_t));
                        fused = true;
                        n_generic_lines++;
                        R_geno_materialised = false;
                        R_geno_line = p;
                        R_geno_end = le;
                        return;
                    }
                    // the pre-parser found a short line: fall through to the general path (lazy retention)
                }
            }
            if (geno_cols && all_included && GT_idx == 0 && fused_walk(p, le, hap, N, scratch.data(), W, sm)) {
                fused = true;
                n_generic_lines++;
                R_geno_materialised = false;
                R_geno_line = p;
                R_geno_end = le;
                return;
            }
#endif
            n_generic_lines++;
            // generic path: walk sample columns of this line
            const char* q = geno_cols ? p : le;
            size_t i = 0;
            if (geno_cols) {
#if SWIFTLD_SIMD
                // tab positions from a 64-byte compare bitmap (masked loads never touch bytes past le)
                const __m512i vtab = _mm512_set1_epi8('\t');
                const char* blk = q;
                auto load_mask = [&](const char* b0) -> uint64_t {
                    const size_t rem = size_t(le - b0);
                    const __mmask64 lm = rem >= 64 ? ~__mmask64(0) : ((__mmask64(1) << rem) - 1);
                    return _mm512_cmpeq_epi8_mask(_mm512_maskz_loadu_epi8(lm, b0), vtab);
                };
                uint64_t tm = load_mask(blk);
                for (; i < N; i++) {
                    while (tm == 0) {
                        blk += 64;
                        if (blk >= le) break;
                        tm = load_mask(blk);
                    }
                    const char* t = tm ? blk + __builtin_ctzll(tm) : nullptr;
                    tm &= tm - 1;
                    const char* fend = t ? t : le;
                    // GT first and exactly 3 characters "a/b" or "a|b" (parse_gt's first branch), else generic
                    if (GT_idx == 0 && fend - q >= 3 && (q[1] == '/' || q[1] == '|') && q[0] != ':' && q[2] != ':' &&
                        (q + 3 == fend || q[3] == ':')) {
                        GT& o = gts[i];
                        o.ploidy = 2;
                        o.phase = q[1];
                        o.a1 = (q[0] != '.') ? (q[0] - '0') : -1;
                        o.a2 = (q[2] != '.') ? (q[2] - '0') : -1;
                        o.poly = false;
                    } else {
                        parse_sample(q, fend, GT_idx, gts[i]);
                    }
                    if (!t) { i++; break; }
                    q = t + 1;
                }
#else
                for (; i < N; i++) {
                    const char* t = (const char*)memchr(q, '\t', size_t(le - q));
                    const char* fend = t ? t : le;
                    parse_sample(q, fend, GT_idx, gts[i]);
                    if (!t) { i++; break; }
                    q = t + 1;
                }
#endif
            }
            if (i < N) {
                // short line: remaining columns keep the previous fully parsed values
                if (!R_geno_materialised) {
                    const char* r = R_geno_line;
                    for (size_t j = 0; j < N; j++) {
                        const char* t = (const char*)memchr(r, '\t', size_t(R_geno_end - r));
                        const char* fend = t ? t : R_geno_end;
                        R_geno[j] = sv(r, size_t(fend - r));
                        r = t ? t + 1 : R_geno_end;
                    }
                    R_geno_materialised = true;
                    R_geno_min = R_geno_line;   // entries point into that line or later ones
                }
                if (geno_cols) {
                    const char* r = p;
                    for (size_t j = 0; j < i; j++) {
                        const char* t = (const char*)memchr(r, '\t', size_t(le - r));
                        const char* fend = t ? t : le;
                        R_geno[j] = sv(r, size_t(fend - r));
                        r = t ? t + 1 : le;
                    }
                }
                for (size_t j = i; j < N; j++) parse_sample(R_geno[j].data(), R_geno[j].data() + R_geno[j].size(), GT_idx, gts[j]);
            } else {
                R_geno_materialised = false;
                R_geno_line = p;
                R_geno_end = le;
            }
            sm = Summ();
            {
                // branch-free (genotypes are data-dependent and mispredict): same sums as
                // entry::get_allele_counts (only alleles 0 and 1 matter for a biallelic site)
                unsigned pa = 0, pi = 0, ua = 0, ui = 0, nnm = 0, nchr = 0;
                long c0 = 0, c1 = 0;
                for (size_t j = 0; j < N; j++) {
                    const GT& g = gts[j];
                    const unsigned inc = unsigned(include[j] != 0), im = 0u - inc;   // im: all ones if included
                    const unsigned po = unsigned(g.poly), un = unsigned(g.phase != '|');
                    pa |= po; pi |= po & inc; ua |= un; ui |= un & inc;
                    c0 += long(im & (unsigned(g.a1 == 0) + unsigned(g.a2 == 0)));
                    c1 += long(im & (unsigned(g.a1 == 1) + unsigned(g.a2 == 1)));
                    nnm += im & (unsigned(g.a1 > -1) + unsigned(g.a2 > -1));
                    nchr += im & unsigned(g.ploidy);
                }
                sm.poly_all = pa; sm.poly_incl = pi; sm.unph_all = ua; sm.unph_incl = ui;
                sm.c0 = c0; sm.c1 = c1; sm.Nnm = nnm; sm.Nchr = nchr;
            }
        };

        // 5. phased_only (filter_sites_by_phase): parses included individuals only
        if (P.phased_only) {
            ensure_parsed();
            if (sm.poly_incl) die("Polyploidy found, and not supported by vcftools");
            if (sm.unph_incl) continue;
        }
        // 6. --maf / --max-missing (filter_sites_by_frequency_and_call_rate): parses all individuals
        if (freq_filter) {
            ensure_parsed();
            if (sm.poly_all) die("Polyploidy found, and not supported by vcftools");
            if (GT_idx == -1) die("Require Genotypes in variant file to filter by frequency and/or call rate");
            double maf = DBL_MAX;
            long cnt[2] = {sm.c0, sm.c1};
            for (int u = 0; u < 2; u++) {
                double freq = int(cnt[u]) / (double)sm.Nnm;
                double folded = std::min(freq, 1.0 - freq);
                maf = std::min(maf, folded);
            }
            if ((maf < P.min_maf) || (maf > P.max_maf)) continue;
            double call_rate = sm.Nnm / double(sm.Nchr);
            if (call_rate < P.min_site_call_rate) continue;
        }
        // 7. --thin (filter_sites_by_thinning)
        if (P.thin >= 1) {
            bool pass = true;
            if (CHROM == sv(thin_chrom)) {
                int dist = int(unsigned(POS) - unsigned(thin_pos));
                if (dist < P.thin) pass = false;
            }
            if (pass) thin_pos = POS;
            thin_chrom = string(CHROM);
            if (!pass) continue;
        }
        // ---- kept site
        N_kept_entries++;
        ensure_parsed();
        if (sm.poly_all) die("Polyploidy found, and not supported by vcftools");
        if (hap && sm.unph_all) die("Require phased haplotypes for r^2 calculation (use --phased)");
        // chrom id
        int cid;
        {
            auto it = chrom_ids.find(string(CHROM));
            if (it == chrom_ids.end()) {
                cid = int(chrom_names.size());
                chrom_ids.emplace(string(CHROM), cid);
                chrom_names.push_back(string(CHROM));
                by_chrom.emplace_back();
                chrom_sorted.push_back(1);
                if (clists.size() >= (1u << 20)) die("too many chromosomes");
                clists.push_back(new ChromList());
            } else cid = it->second;
        }
        const size_t kidx = s_pos.size();
        if (!by_chrom[size_t(cid)].empty() && POS < s_pos[by_chrom[size_t(cid)].back()]) {
            chrom_sorted[size_t(cid)] = 0;
            clists[size_t(cid)]->unsorted.store(true, std::memory_order_release);   // before this site is pushed
        }
        rank_in_chrom.push_back(uint32_t(by_chrom[size_t(cid)].size()));
        by_chrom[size_t(cid)].push_back(uint32_t(kidx));
        s_chrom.push_back(cid);
        s_pos.push_back(POS);
        SiteRec& rec = SS.append();
        uint64_t* pl = SS.planes_mut(kidx);
        rec.pos = POS;
        rec.chrom = uint32_t(cid);
        rec.rank = rank_in_chrom.back();
        {
            char b[24];
            char* e2 = swiftfmt::fmt_i32(b, POS);
            *e2++ = '\t';
            const size_t bl = size_t(e2 - b);
            memcpy(rec.pos2, b, bl);
            rec.pos2len = uint8_t(bl);
            {
                Sites::Block* blk = SS.blocks[kidx / Sites::B];
                char* t = blk->p2n[kidx % Sites::B];
                memcpy(t, b, bl);
                memcpy(t + bl, Hstr, Hlen);   // Hlen <= 11 (digits + tab), bl <= 12
                blk->p2nlen[kidx % Sites::B] = uint8_t(bl + Hlen);
            }
            const size_t plen = CHROM.size() + 1 + bl;
            rec.prelen = uint8_t(plen <= 255 ? plen : 0);
            if (plen <= sizeof(rec.pre)) {
                memcpy(rec.pre, CHROM.data(), CHROM.size());
                rec.pre[CHROM.size()] = '\t';
                memcpy(rec.pre + CHROM.size() + 1, b, bl);
                rec.prelen_long = 0;
            } else {
                string lp = string(CHROM) + "\t" + string(b, bl);
                SS.set_long_prefix(kidx, lp);
                rec.prelen_long = 1;
            }
        }
        // planes
#if SWIFTLD_SIMD
        if (fast) {
            fast_geno_planes(p, size_t(le - p), hap, pl, W);
        } else if (fused) {
            memcpy(pl, scratch.data(), nplanes * W * sizeof(uint64_t));
        } else
#endif
        if (hap) {
            // words built in registers, branch-free (same bits as setting them one by one)
            uint64_t* Pp = pl; uint64_t* Ap = pl + W;
            for (size_t w = 0; w < W; w++) {
                uint64_t P = 0, A = 0;
                const size_t j0 = 32 * w, j1 = std::min(NI, j0 + 32);
                for (size_t j = j0; j < j1; j++) {
                    const GT& g = gts[size_t(incl_idx[j])];
                    int d1, d2;
                    enc_dec(g.a1, g.a2, d1, d2);   // ploidy <= 2 here
                    const unsigned s0 = unsigned(2 * (j - j0));
                    P |= (uint64_t(d1 >= 0) << s0) | (uint64_t(d2 >= 0) << (s0 + 1));
                    A |= (uint64_t(d1 > 0) << s0) | (uint64_t(d2 > 0) << (s0 + 1));
                }
                Pp[w] |= P; Ap[w] |= A;
            }
        } else {
            uint64_t* Vp = pl; uint64_t* ap = pl + W; uint64_t* bp = pl + 2 * W;
            for (size_t w = 0; w < W; w++) {
                uint64_t V = 0, Am = 0, Bm = 0;
                const size_t j0 = 64 * w, j1 = std::min(NI, j0 + 64);
                for (size_t j = j0; j < j1; j++) {
                    const GT& g = gts[size_t(incl_idx[j])];
                    int d1, d2;
                    enc_dec(g.a1, g.a2, d1, d2);
                    // "LD: Only using diploid individuals." (ploidy 2) and both decoded alleles present
                    const uint64_t v = uint64_t((g.ploidy == 2) & (d1 >= 0) & (d2 >= 0));
                    const uint64_t x1 = uint64_t((d1 != d2) | (d1 == 0));   // x = (d1==d2) ? (d1==0 ? 2 : 0) : 1; x >= 1
                    const uint64_t x2 = uint64_t((d1 == d2) & (d1 == 0));   // x == 2
                    const unsigned sh = unsigned(j - j0);
                    V |= v << sh; Am |= (v & x1) << sh; Bm |= (v & x2) << sh;
                }
                Vp[w] |= V; ap[w] |= Am; bp[w] |= Bm;
            }
        }
        if (hap) {
            uint32_t present = popcnt_words(pl, W), alt = popcnt_words(pl + W, W);
            rec.full = present == H;
            rec.c1 = alt;
            rec.c2 = 0;
        } else {
            uint32_t valid = popcnt_words(pl, W), sa = popcnt_words(pl + W, W), sb = popcnt_words(pl + 2 * W, W);
            rec.full = valid == NI;
            rec.c1 = sa + sb;
            rec.c2 = sa + 3 * sb;
        }
        SS.count = kidx + 1;
        clists[size_t(cid)]->push(uint32_t(kidx));
        advance_rows(false);
    }
    if (N_kept_entries <= 1) die("Insufficient sites remained after filtering");
    advance_rows(true);
    if (cur) { Q.publish(cur); cur = nullptr; }
    Q.parse_done.store(true, std::memory_order_release);
    double t_parse = now_s();
    flush_ready();
    worker_loop();
    for (auto& t : extra) t.join();
    if (in.st_thread) in.st_thread->join();
    // drain (a chunk finished while another thread held the writer role, or writes were gated)
    flush_ready();
    if (!header_written) { write_all(out_fd, hdr, strlen(hdr)); header_written = true; out_bytes += strlen(hdr); }
    if (next_write < Q.published.load()) die("internal: unwritten chunks");
    if (!P.stream_out) {
        if (falloc_end != SIZE_MAX && falloc_end > out_bytes) { if (ftruncate(out_fd, off_t(out_bytes))) {} }
        close(out_fd);
    }
    if (prof) {
        double t_end = now_s();
        if (in.st) fprintf(stderr, "[profile] background inflate %.3f s, %zu chunks of 2 MB recycled\n", in.st->t_inflate, in.st->n_recycled.load());
        fprintf(stderr, "[profile] read %.3f s  parse+filter(+overlapped pairs) %.3f s (%zu sites kept, %zu chunks)  rest %.3f s  total %.3f s  [genotype lines: %zu SIMD, %zu generic]  thread-sum format %.3f s, write %.3f s\n",
                t_read - t_start, t_parse - t_read, s_pos.size(), Q.published.load(), t_end - t_parse, t_end - t_start,
                n_fast_lines, n_generic_lines, ns_fmt.load() * 1e-9, ns_write.load() * 1e-9);
    }
    return 0;
}
