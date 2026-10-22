// fmt.h - exact, fast replacements for `ostream << double` (printf "%.6g") and `ostream << int`.
//
// libstdc++ formats a double with the stream's default flags as vsnprintf("%.*g", 6, v) in the "C" locale.
// glibc rounds the exact binary value half-to-even. This formatter reproduces that byte for byte:
//   * nan/-nan/inf/-inf by sign bit, "0"/"-0" for zeros;
//   * for 1e-22 <= |v| < 1e6 the 6 significant digits are computed EXACTLY with 128-bit integer arithmetic
//     (v = m*2^q, scaled = m*5^s*2^(q+s), integer part + exact remainder -> round half to even);
//   * anything else (huge, tiny, subnormal) falls back to snprintf("%.6g"), which is the definition.
#pragma once
#if defined(__AVX512F__)
#include <immintrin.h>
#endif
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>

namespace swiftfmt {

static const char DIGITS2[201] =
    "00010203040506070809101112131415161718192021222324252627282930313233343536373839"
    "40414243444546474849505152535455565758596061626364656667686970717273747576777879"
    "8081828384858687888990919293949596979899";

static inline uint64_t pow5(int s) {
    static const uint64_t T[28] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull,
        1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull,
        152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull,
        476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull,
        298023223876953125ull, 1490116119384765625ull, 7450580596923828125ull};
    return T[s];
}

// Unsigned decimal.
static inline char* fmt_u32(char* p, uint32_t v) {
    char tmp[12];
    char* e = tmp + 12;
    char* q = e;
    while (v >= 100) {
        uint32_t r = v % 100;
        v /= 100;
        q -= 2;
        memcpy(q, DIGITS2 + 2 * r, 2);
    }
    if (v >= 10) {
        q -= 2;
        memcpy(q, DIGITS2 + 2 * v, 2);
    } else {
        *--q = char('0' + v);
    }
    size_t n = size_t(e - q);
    memcpy(p, q, n);
    return p + n;
}

static inline char* fmt_i32(char* p, int32_t v) {
    if (v < 0) {
        *p++ = '-';
        return fmt_u32(p, uint32_t(0) - uint32_t(v));
    }
    return fmt_u32(p, uint32_t(v));
}

static char* fmt_g_slow(char* p, double v) {
    char tmp[64];
    int n = snprintf(tmp, sizeof tmp, "%g", v);
    memcpy(p, tmp, size_t(n));
    return p + n;
}

// Exact "%.6g" (reference implementation, 128-bit integer rounding).
static inline char* fmt_g_exact(char* p, double v) {
    uint64_t bits;
    memcpy(&bits, &v, 8);
    const uint32_t eb = uint32_t(bits >> 52) & 0x7ff;
    const uint64_t frac = bits & ((1ull << 52) - 1);
    const bool neg = (bits >> 63) != 0;
    if (eb == 0x7ff) {
        if (frac) {
            if (neg) { memcpy(p, "-nan", 4); return p + 4; }
            memcpy(p, "nan", 3); return p + 3;
        }
        if (neg) { memcpy(p, "-inf", 4); return p + 4; }
        memcpy(p, "inf", 3); return p + 3;
    }
    if (eb == 0) {
        if (frac == 0) {
            if (neg) { memcpy(p, "-0", 2); return p + 2; }
            *p = '0'; return p + 1;
        }
        return fmt_g_slow(p, v);  // subnormal
    }
    // v = m * 2^q
    const uint64_t m = frac | (1ull << 52);
    const int q = int(eb) - 1075;
    // floor(log10(2^(eb-1023))) estimate; true E is this or this+1.
    int E = ((int(eb) - 1023) * 78913) >> 18;
    if (E >= 6 || E < -22) return fmt_g_slow(p, v);
    uint32_t R;
    for (;;) {
        const int s = 5 - E;          // 0 <= s <= 27
        const int shift = -(q + s);   // > 0 in the supported range
        const unsigned __int128 N = (unsigned __int128)m * pow5(s);
        const unsigned __int128 R0 = N >> shift;
        if (R0 >= 1000000) { E++; if (E >= 6) return fmt_g_slow(p, v); continue; }
        if (R0 < 100000) { E--; if (E < -22) return fmt_g_slow(p, v); continue; }
        const unsigned __int128 rem = N - (R0 << shift);
        const unsigned __int128 half = (unsigned __int128)1 << (shift - 1);
        R = uint32_t(R0);
        if (rem > half || (rem == half && (R & 1))) R++;
        if (R == 1000000) { R = 100000; E++; }
        break;
    }
    if (neg) *p++ = '-';
    char d[6];
    {
        uint32_t hi = R / 10000, mid = (R / 100) % 100, lo = R % 100;
        memcpy(d, DIGITS2 + 2 * hi, 2);
        memcpy(d + 2, DIGITS2 + 2 * mid, 2);
        memcpy(d + 4, DIGITS2 + 2 * lo, 2);
    }
    int nd = 6;  // significant digits after stripping trailing zeros
    while (nd > 1 && d[nd - 1] == '0') nd--;
    if (E >= 6 || E < -4) {
        *p++ = d[0];
        if (nd > 1) {
            *p++ = '.';
            memcpy(p, d + 1, size_t(nd - 1));
            p += nd - 1;
        }
        *p++ = 'e';
        int ae = E;
        if (ae < 0) { *p++ = '-'; ae = -ae; } else { *p++ = '+'; }
        if (ae < 10) { *p++ = '0'; *p++ = char('0' + ae); }
        else if (ae < 100) { memcpy(p, DIGITS2 + 2 * ae, 2); p += 2; }
        else { p = fmt_u32(p, uint32_t(ae)); }
        return p;
    }
    if (E >= 0) {
        // integer part: d[0..E]
        memcpy(p, d, size_t(E + 1));
        p += E + 1;
        if (nd > E + 1) {
            *p++ = '.';
            memcpy(p, d + E + 1, size_t(nd - (E + 1)));
            p += nd - (E + 1);
        }
        return p;
    }
    // -4 <= E <= -1: "0." + (-E-1) zeros + digits
    *p++ = '0';
    *p++ = '.';
    for (int z = 0; z < -E - 1; z++) *p++ = '0';
    memcpy(p, d, size_t(nd));
    return p + nd;
}


// ---------------------------------------------------------------- fast exact "%.6g"
// E = floor(log10 v) from the binary exponent plus one threshold compare; 6 significant digits from one
// double multiply by an exact power of ten (error <= 1.2e-10 absolute on a value < 1e6). If the fractional part is
// within 1e-6 of a rounding tie, or the estimate is off by one near a power of ten, defer to fmt_g_exact.
struct G6Tables {
    signed char E[2048];
    double thr[2048];
    double pow10[23];
    G6Tables() {
        double p = 1.0;
        for (int s = 0; s <= 22; s++) { pow10[s] = p; p *= 10.0; }
        // 10^k as long double, built by the same repeated *10 / /10 chains as a per-entry loop would (identical
        // values), once for k in [-99, 101] instead of once per exponent (the per-entry loop cost ~0.5 ms at startup)
        long double tk[201];
        tk[99] = 1.0L;
        for (int k = 1; k <= 101; k++) tk[99 + k] = tk[99 + k - 1] * 10.0L;
        for (int k = 1; k <= 99; k++) tk[99 - k] = tk[99 - k + 1] / 10.0L;
        for (int eb = 0; eb < 2048; eb++) {
            int e2 = eb - 1023;
            int E0 = (e2 * 78913) >> 18;  // floor(e2 * log10(2)), exact for |e2| < 1650
            if (E0 < -100) E0 = -100;
            if (E0 > 100) E0 = 100;
            E[eb] = (signed char)E0;
            // 10^(E0+1) as the nearest double; +inf if it is not inside this binade
            const long double t = tk[99 + E0 + 1];
            long double hi = ldexpl(1.0L, e2 + 1);
            thr[eb] = (t < hi) ? double(t) : __builtin_inf();
        }
    }
};
static const G6Tables G6;

// 3-digit ASCII groups, little-endian in a uint32 (byte 0 = hundreds digit)
struct Lut3 {
    uint32_t t[1000];
    Lut3() {
        for (uint32_t i = 0; i < 1000; i++)
            t[i] = uint32_t('0' + i / 100) | (uint32_t('0' + (i / 10) % 10) << 8) | (uint32_t('0' + i % 10) << 16);
    }
};
static const Lut3 L3;

static inline char* fmt_g(char* p, double v) {
    uint64_t bits;
    memcpy(&bits, &v, 8);
    const uint32_t eb = uint32_t(bits >> 52) & 0x7ff;
    if (__builtin_expect(eb == 0 || eb == 0x7ff, 0)) return fmt_g_exact(p, v);
    const double a = __builtin_fabs(v);
    int E = G6.E[eb] + (a >= G6.thr[eb]);
    if (__builtin_expect(E < -17 || E > 5, 0)) return fmt_g_exact(p, v);
    const double y = a * G6.pow10[5 - E];
    uint32_t R = uint32_t(y);
    const double fr = y - double(R);
    if (__builtin_expect(R < 100000 || R >= 1000000 || __builtin_fabs(fr - 0.5) < 1e-6, 0)) return fmt_g_exact(p, v);
    R += (fr > 0.5);
    *p = '-';
    p += (bits >> 63);
    if (__builtin_expect(R == 1000000, 0)) { R = 100000; E++; }
    const uint32_t qh = R / 1000, ql = R - qh * 1000;
    const uint64_t dg = uint64_t(L3.t[qh]) | (uint64_t(L3.t[ql]) << 24);   // 6 ASCII digits, bytes 6-7 zero
    const uint64_t x = dg ^ 0x0000303030303030ull;                          // nonzero bytes = digits != '0'
    const int nd = ((63 - __builtin_clzll(x)) >> 3) + 1;                    // digits left after stripping zeros
    if (E < 0 && E >= -4) {
        const int z = -E - 1;
        memcpy(p, "0.000000", 8);
        memcpy(p + 2 + z, &dg, 8);
        return p + 2 + z + nd;
    }
    const uint64_t rest = dg >> 8;   // digits 1..5
    if (E == 0) {
        p[0] = char(dg);
        p[1] = '.';
        memcpy(p + 2, &rest, 8);
        return p + (nd > 1 ? nd + 1 : 1);
    }
    if (E < 0 || E >= 6) {
        p[0] = char(dg);
        p[1] = '.';
        memcpy(p + 2, &rest, 8);
        p += (nd > 1 ? nd + 1 : 1);
        int ae = E < 0 ? -E : E;
        p[0] = 'e';
        p[1] = E < 0 ? '-' : '+';
        memcpy(p + 2, DIGITS2 + 2 * ae, 2);  // |E| <= 17 here (E == 6 only after a carry)
        return p + 4;
    }
    // 1 <= E <= 5
    char d[8];
    memcpy(d, &dg, 8);
    memcpy(p, d, size_t(E + 1));
    p += E + 1;
    if (nd > E + 1) {
        *p++ = '.';
        memcpy(p, d + E + 1, size_t(nd - (E + 1)));
        p += nd - (E + 1);
    }
    return p;
}


// ---------------------------------------------------------------- batch digit generation (AVX-512)
// g6_batch converts n doubles to (6 ASCII digits, exponent, digit count, sign, ok) with 8 lanes at a time.
// Lanes that are not "ok" (zero, subnormal, inf/nan, |v| outside [1e-11, 1e5), or within 1e-6 of a rounding
// tie) must be printed with fmt_g, which is exact. emit_g6 prints one value from its batch result.
#if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512CD__) && defined(__AVX512VL__)
#define SWIFTFMT_BATCH 1
struct G6Batch {
    uint64_t dg;     // 6 ASCII digits, byte 0 = most significant
    uint32_t meta;   // byte0: E+64, byte1: nd, byte2: ok, byte3: negative
};

static inline void g6_batch(const double* v, size_t n, uint64_t* dgo, uint32_t* mo) {
    alignas(64) static const double P10[16] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                               1e12, 1e13, 1e14, 1e15};
    const __m512d plo = _mm512_load_pd(P10), phi = _mm512_load_pd(P10 + 8);
    const __m512i c7ff = _mm512_set1_epi64(0x7ff), c1023 = _mm512_set1_epi64(1023);
    const __m512i c78913 = _mm512_set1_epi64(78913);
    const __m512d absmask = _mm512_castsi512_pd(_mm512_set1_epi64(0x7fffffffffffffffll));
    const __m512d c1e6 = _mm512_set1_pd(1e6), chalf = _mm512_set1_pd(0.5), ctie = _mm512_set1_pd(1e-6);
    const __m512i c5 = _mm512_set1_epi64(5), c0 = _mm512_setzero_si512(), c15 = _mm512_set1_epi64(15);
    const __m512i c1 = _mm512_set1_epi64(1);
    const __m512i r_lo = _mm512_set1_epi64(100000), r_hi = _mm512_set1_epi64(1000000);
    const __m512i k1000 = _mm512_set1_epi64(1000), k4294968 = _mm512_set1_epi64(4294968);
    const __m512i k41 = _mm512_set1_epi64(41), k100 = _mm512_set1_epi64(100), k103 = _mm512_set1_epi64(103);
    const __m512i k10 = _mm512_set1_epi64(10), kzero = _mm512_set1_epi64(0x303030303030ll);
    const __m512i k64 = _mm512_set1_epi64(64);
    for (size_t i = 0; i < n; i += 8) {
        const __mmask8 lm = (n - i >= 8) ? 0xff : __mmask8((1u << (n - i)) - 1);
        const __m512d x = _mm512_maskz_loadu_pd(lm, v + i);
        const __m512i bits = _mm512_castpd_si512(x);
        const __m512i eb = _mm512_and_si512(_mm512_srli_epi64(bits, 52), c7ff);
        __mmask8 ok = lm & _mm512_cmpneq_epi64_mask(eb, c0) & _mm512_cmpneq_epi64_mask(eb, c7ff);
        const __m512d a = _mm512_and_pd(x, absmask);
        const __m512i e2 = _mm512_sub_epi64(eb, c1023);
        const __m512i E0 = _mm512_srai_epi64(_mm512_mul_epi32(e2, c78913), 18);
        __m512i s0 = _mm512_sub_epi64(c5, E0);                       // want 1..15
        ok &= _mm512_cmpge_epi64_mask(s0, c1) & _mm512_cmple_epi64_mask(s0, c15);
        s0 = _mm512_min_epi64(_mm512_max_epi64(s0, c1), c15);
        const __m512d y0 = _mm512_mul_pd(a, _mm512_permutex2var_pd(plo, s0, phi));
        const __m512d y1 = _mm512_mul_pd(a, _mm512_permutex2var_pd(plo, _mm512_sub_epi64(s0, c1), phi));
        const __mmask8 big = _mm512_cmp_pd_mask(y0, c1e6, _CMP_GE_OQ);
        const __m512d y = _mm512_mask_blend_pd(big, y0, y1);
        __m512i E = _mm512_mask_add_epi64(E0, big, E0, c1);
        __m512i R = _mm512_cvttpd_epu64(y);
        const __m512d fr = _mm512_sub_pd(y, _mm512_cvtepu64_pd(R));
        ok &= _mm512_cmpge_epu64_mask(R, r_lo) & _mm512_cmplt_epu64_mask(R, r_hi);
        ok &= _mm512_cmp_pd_mask(_mm512_and_pd(_mm512_sub_pd(fr, chalf), absmask), ctie, _CMP_GE_OQ);
        R = _mm512_mask_add_epi64(R, _mm512_cmp_pd_mask(fr, chalf, _CMP_GT_OQ), R, c1);
        const __mmask8 carry = _mm512_cmpeq_epi64_mask(R, r_hi);
        R = _mm512_mask_mov_epi64(R, carry, r_lo);
        E = _mm512_mask_add_epi64(E, carry, E, c1);
        // digits
        const __m512i qh = _mm512_srli_epi64(_mm512_mul_epu32(R, k4294968), 32);
        const __m512i ql = _mm512_sub_epi64(R, _mm512_mul_epu32(qh, k1000));
        const __m512i h0 = _mm512_srli_epi64(_mm512_mul_epu32(qh, k41), 12);
        const __m512i hr = _mm512_sub_epi64(qh, _mm512_mul_epu32(h0, k100));
        const __m512i h1 = _mm512_srli_epi64(_mm512_mul_epu32(hr, k103), 10);
        const __m512i h2 = _mm512_sub_epi64(hr, _mm512_mul_epu32(h1, k10));
        const __m512i l0 = _mm512_srli_epi64(_mm512_mul_epu32(ql, k41), 12);
        const __m512i lr = _mm512_sub_epi64(ql, _mm512_mul_epu32(l0, k100));
        const __m512i l1 = _mm512_srli_epi64(_mm512_mul_epu32(lr, k103), 10);
        const __m512i l2 = _mm512_sub_epi64(lr, _mm512_mul_epu32(l1, k10));
        __m512i dg = _mm512_or_si512(h0, _mm512_slli_epi64(h1, 8));
        dg = _mm512_or_si512(dg, _mm512_slli_epi64(h2, 16));
        dg = _mm512_or_si512(dg, _mm512_slli_epi64(l0, 24));
        dg = _mm512_or_si512(dg, _mm512_slli_epi64(l1, 32));
        dg = _mm512_or_si512(dg, _mm512_slli_epi64(l2, 40));
        dg = _mm512_add_epi64(dg, kzero);
        const __m512i nz = _mm512_xor_si512(dg, kzero);
        const __m512i nd = _mm512_add_epi64(_mm512_srli_epi64(_mm512_sub_epi64(_mm512_set1_epi64(63),
                                                                               _mm512_lzcnt_epi64(nz)), 3), c1);
        const __m512i neg = _mm512_srli_epi64(bits, 63);
        __m512i meta = _mm512_and_si512(_mm512_add_epi64(E, k64), _mm512_set1_epi64(0xff));
        meta = _mm512_or_si512(meta, _mm512_slli_epi64(nd, 8));
        meta = _mm512_or_si512(meta, _mm512_maskz_mov_epi64(ok, _mm512_set1_epi64(1 << 16)));
        meta = _mm512_or_si512(meta, _mm512_slli_epi64(neg, 24));
        _mm512_mask_storeu_epi64(dgo + i, lm, dg);
        _mm256_mask_storeu_epi32(mo + i, lm, _mm512_cvtepi64_epi32(meta));
    }
}

// Branch-free layout: the 6 digits plus the constant characters are placed in one 16-byte vector
//   bytes 0-5 digits, 6 '0', 7 '.', 8 'e', 9 exponent sign, 10-11 exponent digits
// and a per-(E, nd) shuffle mask produces the %g text; a length table gives how far to advance.
struct G6Layout {
    alignas(16) uint8_t mask[32][8][16];   // [E + 20][nd]
    uint8_t len[32][8];
    uint16_t expd[32];                     // two exponent digits for |E|
    G6Layout() {
        memset(mask, 0x80, sizeof mask);
        memset(len, 0, sizeof len);
        for (int E = -20; E <= 11; E++) {
            int ae = E < 0 ? -E : E;
            expd[E + 20] = uint16_t(uint8_t('0' + (ae / 10) % 10) | (uint8_t('0' + ae % 10) << 8));
            for (int nd = 1; nd <= 6; nd++) {
                uint8_t m[16];
                int k = 0;
                if (E >= -4 && E < 0) {
                    m[k++] = 6; m[k++] = 7;
                    for (int z = 0; z < -E - 1; z++) m[k++] = 6;
                    for (int d = 0; d < nd; d++) m[k++] = uint8_t(d);
                } else if (E >= 0 && E <= 5) {
                    for (int d = 0; d <= E; d++) m[k++] = uint8_t(d);
                    if (nd > E + 1) { m[k++] = 7; for (int d = E + 1; d < nd; d++) m[k++] = uint8_t(d); }
                } else {
                    m[k++] = 0;
                    if (nd > 1) { m[k++] = 7; for (int d = 1; d < nd; d++) m[k++] = uint8_t(d); }
                    m[k++] = 8; m[k++] = 9; m[k++] = 10; m[k++] = 11;
                }
                memcpy(mask[E + 20][nd], m, size_t(k));
                len[E + 20][nd] = uint8_t(k);
            }
        }
    }
};
static const G6Layout GL;

static inline char* emit_g6(char* p, double v, uint64_t dg, uint32_t meta) {
    if (__builtin_expect(!(meta & (1u << 16)), 0)) return fmt_g(p, v);
    *p = '-';
    p += (meta >> 24);
    const int ei = int(meta & 0xff) - 64 + 20;   // E + 20, in [10, 26] on the fast path
    const int nd = int((meta >> 8) & 0xff);
    const uint64_t lo = dg | (0x30ull << 48) | (0x2eull << 56);
    const uint64_t hiw = uint64_t('e') | (uint64_t(ei < 20 ? '-' : '+') << 8) | (uint64_t(GL.expd[ei]) << 16);
    const __m128i src = _mm_set_epi64x((long long)hiw, (long long)lo);
    const __m128i out = _mm_shuffle_epi8(src, _mm_load_si128((const __m128i*)GL.mask[ei][nd]));
    _mm_storeu_si128((__m128i*)p, out);
    return p + GL.len[ei][nd];
}

// ---- variant with the whole field in one 16-byte source: g6_batch16 writes per value
//   bytes 0-5 digits, 6 '0', 7 '.', 8 'e', 9 exponent sign, 10-11 exponent digits, 12 '-', 13 '\t'
// and a 16-bit code: (neg << 8) | ((E + 20) << 3) | nd, or 0x8000 when the lane must take the exact path.
// emit16 prints the value followed by '\t' with one shuffle (the sign and the tab are selected by the mask).
struct G6Layout16 {
    alignas(16) uint8_t mask[512][16];
    uint8_t len[512];
    G6Layout16() {
        memset(mask, 0x80, sizeof mask);
        memset(len, 0, sizeof len);
        // special codes (nd = 0, unused by finite values): ei 0 = zero ("0"), 1 = NaN ("nan"), 2 = inf ("inf");
        // src bytes 0-2 hold the text, 12 '-', 13 '\t'
        for (int neg = 0; neg < 2; neg++)
            for (int sp = 0; sp < 3; sp++) {
                const int code = (neg << 8) | (sp << 3);
                uint8_t* m = mask[code];
                int k = 0;
                if (neg) m[k++] = 12;
                const int nch = sp == 0 ? 1 : 3;
                for (int t = 0; t < nch; t++) m[k++] = uint8_t(t);
                m[k++] = 13;
                len[code] = uint8_t(k);
            }
        for (int neg = 0; neg < 2; neg++)
            for (int ei = 0; ei < 32; ei++)
                for (int nd = 1; nd <= 6; nd++) {
                    const int code = (neg << 8) | (ei << 3) | nd;
                    uint8_t* m = mask[code];
                    int k = 0;
                    if (neg) m[k++] = 12;
                    for (int t = 0; t < GL.len[ei][nd]; t++) m[k++] = GL.mask[ei][nd][t];
                    m[k++] = 13;
                    len[code] = uint8_t(k);
                }
    }
};
static const G6Layout16 GL16;

// Exponent text of the 16-byte field source for E + 20 in [0, 32): 'e', sign, two digits of |E|, '-', '\t'
// (the same bytes the arithmetic version produced); looked up with two permutes instead of ~11 vector ops.
struct G6ExpText {
    alignas(64) uint64_t t[32];
    constexpr G6ExpText() : t() {
        for (int i = 0; i < 32; i++) {
            const int E = i - 20, ae = E < 0 ? -E : E;
            t[i] = 0x65ull | (uint64_t(E < 0 ? '-' : '+') << 8) | (uint64_t('0' + ae / 10) << 16) |
                   (uint64_t('0' + ae % 10) << 24) | (0x2dull << 32) | (0x09ull << 40);
        }
    }
};
static constexpr G6ExpText G6EXP{};

static inline void g6_batch16(const double* v, size_t n, uint8_t* srco, uint16_t* codeo) {
    alignas(64) static const double P10[16] = {1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11,
                                               1e12, 1e13, 1e14, 1e15};
    const __m512d plo = _mm512_load_pd(P10), phi = _mm512_load_pd(P10 + 8);
    const __m512i c7ff = _mm512_set1_epi64(0x7ff), c1023 = _mm512_set1_epi64(1023);
    const __m512i c78913 = _mm512_set1_epi64(78913);
    const __m512d absmask = _mm512_castsi512_pd(_mm512_set1_epi64(0x7fffffffffffffffll));
    const __m512d c1e6 = _mm512_set1_pd(1e6), chalf = _mm512_set1_pd(0.5), ctie = _mm512_set1_pd(1e-6);
    const __m512i c5 = _mm512_set1_epi64(5), c4 = _mm512_set1_epi64(4), c0 = _mm512_setzero_si512(), c15 = _mm512_set1_epi64(15);
    const __m512i c1 = _mm512_set1_epi64(1);
    const __m512i r_lo = _mm512_set1_epi64(100000), r_hi = _mm512_set1_epi64(1000000);
    const __m512i k1000 = _mm512_set1_epi64(1000), k4294968 = _mm512_set1_epi64(4294968);
    const __m512i w41 = _mm512_set1_epi16(41), w100 = _mm512_set1_epi16(100), w103 = _mm512_set1_epi16(103),
                  w10 = _mm512_set1_epi16(10);
    const __m512i dshuf = _mm512_broadcast_i32x4(_mm_setr_epi8(0, 1, 4, 2, 3, 6, -128, -128,
                                                               8, 9, 12, 10, 11, 14, -128, -128));
    const __m512i kzero = _mm512_set1_epi64(0x2e30303030303030ll);           // digits '0' + "0."
    const __m512i k20 = _mm512_set1_epi64(20);
    // high word: 'e', exponent sign, two exponent digits, '-', '\t'
    const __m512i khi = _mm512_set1_epi64(0x092d30302d65ll);
    const __m512i xt0 = _mm512_load_si512((const void*)G6EXP.t), xt1 = _mm512_load_si512((const void*)(G6EXP.t + 8)),
                  xt2 = _mm512_load_si512((const void*)(G6EXP.t + 16)), xt3 = _mm512_load_si512((const void*)(G6EXP.t + 24));
    const __m512i c16 = _mm512_set1_epi64(16);
    const __m512i ix_lo = _mm512_set_epi64(11, 3, 10, 2, 9, 1, 8, 0), ix_hi = _mm512_set_epi64(15, 7, 14, 6, 13, 5, 12, 4);
    for (size_t i = 0; i < n; i += 8) {
        const size_t rem = n - i;
        const __mmask8 lm = (rem >= 8) ? 0xff : __mmask8((1u << rem) - 1);
        const __m512d x = _mm512_maskz_loadu_pd(lm, v + i);
        const __m512i bits = _mm512_castpd_si512(x);
        const __m512i eb = _mm512_and_si512(_mm512_srli_epi64(bits, 52), c7ff);
        const __mmask8 e0 = _mm512_cmpeq_epi64_mask(eb, c0), emax = _mm512_cmpeq_epi64_mask(eb, c7ff);
        __mmask8 ok = lm & __mmask8(~(e0 | emax));
        const __m512d a = _mm512_and_pd(x, absmask);
        const __m512i e2 = _mm512_sub_epi64(eb, c1023);
        const __m512i E0 = _mm512_srai_epi64(_mm512_mul_epi32(e2, c78913), 18);
        // scale 10^s0 with s0 = 5 - E0 in [1, 15] (s0 - 1 = 4 - E0 in [0, 14], one unsigned compare); the permutes
        // use index bits 3:0 only, and lanes outside the range are not ok, so their results are never used
        const __m512i s0 = _mm512_sub_epi64(c5, E0), s0m1 = _mm512_sub_epi64(c4, E0);
        ok &= _mm512_cmplt_epu64_mask(s0m1, c15);
        const __m512d y0 = _mm512_mul_pd(a, _mm512_permutex2var_pd(plo, s0, phi));
        const __m512d y1 = _mm512_mul_pd(a, _mm512_permutex2var_pd(plo, s0m1, phi));
        const __mmask8 big = _mm512_cmp_pd_mask(y0, c1e6, _CMP_GE_OQ);
        const __m512d y = _mm512_mask_blend_pd(big, y0, y1);
        __m512i E = _mm512_mask_add_epi64(E0, big, E0, c1);
        __m512i R = _mm512_cvttpd_epu64(y);
        const __m512d fr = _mm512_sub_pd(y, _mm512_cvtepu64_pd(R));
        ok &= _mm512_cmpge_epu64_mask(R, r_lo) & _mm512_cmplt_epu64_mask(R, r_hi);
        ok &= _mm512_cmp_pd_mask(_mm512_and_pd(_mm512_sub_pd(fr, chalf), absmask), ctie, _CMP_GE_OQ);
        R = _mm512_mask_add_epi64(R, _mm512_cmp_pd_mask(fr, chalf, _CMP_GT_OQ), R, c1);
        const __mmask8 carry = _mm512_cmpeq_epi64_mask(R, r_hi);
        R = _mm512_mask_mov_epi64(R, carry, r_lo);
        E = _mm512_mask_add_epi64(E, carry, E, c1);
        const __m512i Ec = _mm512_add_epi64(E, k20);
        ok &= _mm512_cmplt_epu64_mask(Ec, _mm512_set1_epi64(32));   // 0 <= Ec < 32 (unsigned compare)
        // digits: R = hi * 1000 + lo; both 3-digit groups converted at once in 16-bit lanes
        // (d0 = x*41 >> 12, d1 = r*103 >> 10 are exact for x, r < 1000), then one byte shuffle
        const __m512i qh = _mm512_srli_epi64(_mm512_mul_epu32(R, k4294968), 32);       // R / 1000
        const __m512i ql = _mm512_sub_epi64(R, _mm512_mul_epu32(qh, k1000));           // R % 1000
        const __m512i x16 = _mm512_or_si512(qh, _mm512_slli_epi64(ql, 16));            // 16-bit lanes [hi, lo, 0, 0]
        const __m512i d0 = _mm512_srli_epi16(_mm512_mullo_epi16(x16, w41), 12);
        const __m512i r16 = _mm512_sub_epi16(x16, _mm512_mullo_epi16(d0, w100));
        const __m512i d1 = _mm512_srli_epi16(_mm512_mullo_epi16(r16, w103), 10);
        const __m512i d2 = _mm512_sub_epi16(r16, _mm512_mullo_epi16(d1, w10));
        // bytes [d0h d1h d0l d1l d2h 0 d2l 0] -> [d0h d1h d2h d0l d1l d2l 0 0]
        const __m512i w = _mm512_or_si512(_mm512_or_si512(d0, _mm512_slli_epi16(d1, 8)), _mm512_slli_epi64(d2, 32));
        __m512i dg = _mm512_shuffle_epi8(w, dshuf);
        // nd: number of significant digits after stripping trailing zeros (digits are 0..9 before adding '0')
        const __m512i nd = _mm512_add_epi64(_mm512_srli_epi64(_mm512_sub_epi64(_mm512_set1_epi64(63),
                                                                               _mm512_lzcnt_epi64(dg)), 3), c1);
        dg = _mm512_add_epi64(dg, kzero);
        // exponent text
        // (lanes with Ec outside [0, 32) are not ok and never emitted from this record)
        __m512i hi = _mm512_mask_blend_epi64(_mm512_cmpge_epi64_mask(Ec, c16),
                                             _mm512_permutex2var_epi64(xt0, Ec, xt1), _mm512_permutex2var_epi64(xt2, Ec, xt3));
        const __m512i neg = _mm512_srli_epi64(bits, 63);
        __m512i code = _mm512_or_si512(_mm512_slli_epi64(neg, 8), _mm512_slli_epi64(Ec, 3));
        code = _mm512_or_si512(code, nd);
        code = _mm512_mask_mov_epi64(_mm512_set1_epi64(0x8000), ok, code);
        // zeros, NaNs and infinities are printed from fixed texts ("0", "nan", "inf", with '-' from the sign bit)
        if (__builtin_expect((e0 | emax) & lm, 0)) {
            const __m512i mant = _mm512_and_si512(bits, _mm512_set1_epi64(0x000fffffffffffffll));
            const __mmask8 mz = _mm512_cmpeq_epi64_mask(mant, c0);
            const __mmask8 is_zero = lm & e0 & mz, is_nan = lm & emax & ~mz, is_inf = lm & emax & mz;
            const __m512i nsh = _mm512_slli_epi64(neg, 8);
            dg = _mm512_mask_mov_epi64(dg, is_zero, _mm512_set1_epi64('0'));
            dg = _mm512_mask_mov_epi64(dg, is_nan, _mm512_set1_epi64(0x6e616e));   // "nan"
            dg = _mm512_mask_mov_epi64(dg, is_inf, _mm512_set1_epi64(0x666e69));   // "inf"
            code = _mm512_mask_mov_epi64(code, is_zero, nsh);
            code = _mm512_mask_mov_epi64(code, is_nan, _mm512_or_si512(nsh, _mm512_set1_epi64(1 << 3)));
            code = _mm512_mask_mov_epi64(code, is_inf, _mm512_or_si512(nsh, _mm512_set1_epi64(2 << 3)));
            hi = _mm512_mask_mov_epi64(hi, is_zero | is_nan | is_inf, khi);   // the exponent bytes are garbage there
        }
        // interleave (lo, hi) into 16-byte records
        const __m512i r0 = _mm512_permutex2var_epi64(dg, ix_lo, hi);   // records 0-3
        const __m512i r1 = _mm512_permutex2var_epi64(dg, ix_hi, hi);   // records 4-7
        uint8_t* so = srco + 16 * i;
        if (rem >= 8) {
            _mm512_storeu_si512(so, r0);
            _mm512_storeu_si512(so + 64, r1);
        } else {
            const uint32_t q = uint32_t(2 * rem);   // qwords to store
            _mm512_mask_storeu_epi64(so, __mmask8(q >= 8 ? 0xff : (1u << q) - 1), r0);
            if (q > 8) _mm512_mask_storeu_epi64(so + 64, __mmask8((1u << (q - 8)) - 1), r1);
        }
        _mm_mask_storeu_epi16(codeo + i, lm, _mm512_cvtepi64_epi16(code));
    }
}

static inline char* emit16(char* p, double v, const uint8_t* src, uint16_t code) {
    if (__builtin_expect(code & 0x8000, 0)) { p = fmt_g(p, v); *p++ = '\t'; return p; }
    const __m128i out = _mm_shuffle_epi8(_mm_loadu_si128((const __m128i*)src),
                                         _mm_load_si128((const __m128i*)GL16.mask[code]));
    _mm_storeu_si128((__m128i*)p, out);
    return p + GL16.len[code];
}

#else
#define SWIFTFMT_BATCH 0
#endif

}  // namespace swiftfmt
