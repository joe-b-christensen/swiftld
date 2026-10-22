# swiftld

swiftld is a drop-in replacement for the pairwise linkage disequilibrium modes of
[vcftools](https://github.com/vcftools/vcftools) 0.1.17, `--hap-r2` and `--geno-r2`. It writes the same `.ld` file,
byte for byte, much faster.

On three 1000 Genomes chromosome 22 regions (503 to 2,504 samples) it was 826× faster than vcftools and 7.6× faster
than PLINK 2, the fastest alternative, which writes a different file (geometric means, 2 threads, AVX-512 CPU; PLINK 2
timed from VCF input to LD output). One megabase with 503 samples took vcftools 29 s and swiftld 0.049 s. Without
AVX-512 the gain is smaller: 102× at 12 threads on an AVX2 laptop. Details are in the preprint (see Citation).

## Build

Requirements: Linux on x86-64, g++ ≥ 9 with C++17, zlib headers and a static zlib (`zlib1g-dev` on Debian/Ubuntu).
Tested with gcc 13.3.0 on Ubuntu 24.04.

```
bash solver/build.sh        # produces solver/swiftld
```

The build uses `-march=native`, so the binary is tuned for the build machine's CPU. The fastest code paths need
AVX-512; on other x86-64 CPUs portable paths give the same output, more slowly.

## Usage

swiftld takes the vcftools options below, in any order, and writes the same output:

```
solver/swiftld --gzvcf in.vcf.gz --hap-r2 --maf 0.05 --ld-window-bp 200000 --out out    # -> out.hap.ld
SWIFTLD_THREADS=4 solver/swiftld --vcf in.vcf --geno-r2 --stdout > out.geno.ld
```

| Group | Options |
|---|---|
| Input / output | `--vcf FILE`, `--gzvcf FILE`, `--out PREFIX`, `--stdout` (also `-c`), `--temp DIR` (accepted, unused) |
| Mode | `--hap-r2` or `--geno-r2` (exactly one) |
| Windows | `--ld-window N`, `--ld-window-min N`, `--ld-window-bp N`, `--ld-window-bp-min N`, `--min-r2 X` |
| Site filters | `--chr C`, `--maf X`, `--max-missing X`, `--remove-filtered-all`, `--remove-indels`, `--thin N`, `--phased` |
| Individual filters | `--keep FILE`, `--remove FILE` |

- The `.ld` output (or standard output with `--stdout`) is byte-identical to vcftools 0.1.17 with the same options,
  at every thread count. Where vcftools fails on an input, swiftld exits non-zero as well.
- Any other option exits with status 2. The `.log` file is not written; messages go to standard error.
- `SWIFTLD_THREADS` sets the number of threads (default: all CPUs available to the process).
- Input: plain or gzip-compressed VCF.

## Check it yourself

```
python3 prepare.py              # builds vcftools 0.1.17 (commit 1c53c3c), generates the synthetic sets and
                                # records vcftools' outputs (needs gcc, g++, zlib, numpy, network; about 25 min)
bash solver/build.sh
python3 evaluate.py             # 300 random VCF files and six benchmark sets compared byte for byte, plus timing
python3 evaluate.py --holdout   # the five holdout sets
```

`data/` holds the five 1000 Genomes regions used in the preprint (see `data/README.txt`).

## Licence

LGPL-3.0 (`LICENSE`), as a derivative of vcftools. Vendored code: libdeflate (MIT); see `LICENSE` and
`solver/LICENSE_NOTE.txt`.

## Citation

Christensen J. SwiftLD: an agent-built, byte-identical vcftools LD that is 826-fold faster than vcftools and 7.6-fold
faster than PLINK 2. Preprint, 2026. Please also cite vcftools (Danecek P et al., Bioinformatics 27:2156–2158,
2011). See `CITATION.cff`.
