#!/bin/bash
# Builds solver/swiftld. evaluate.py runs it from the repository root.
# -ffp-contract=off: no FMA contraction, so every double expression rounds exactly like the reference (-O2, SSE2).
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p .obj
for f in deflate_decompress gzip_decompress crc32 utils; do
  gcc -O2 -march=native -w -c libdeflate/lib/$f.c -o .obj/$f.o &
done
gcc -O2 -march=native -w -c libdeflate/lib/x86/cpu_features.c -o .obj/cpu_features.o &
wait
g++ -std=c++17 -O3 -march=native -mprefer-vector-width=512 -ffp-contract=off -fno-math-errno -Wall -Wno-unused-result \
    -o swiftld swiftld.cpp .obj/*.o -static -lz -pthread
