#!/usr/bin/env bash
# External C baselines on exactly the same inputs: parasail (SIMD NW/SW) and FreeSASA (Shrake-Rupley).
# Builds both into .external/, dumps our inputs and scores, runs them, and checks that parasail's scores are
# identical to ours (BioJava gap open 10 + extend 1 is parasail open 11, extend 1).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
X="$ROOT/.external"; mkdir -p "$X/sets" "$X/asa"
cd "$X"
if [[ ! -f parasail/build/libparasail.a ]]; then
  git clone -q --depth 1 https://github.com/jeffdaily/parasail.git
  (cd parasail && mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF .. >/dev/null && make -j"$(nproc)" parasail >/dev/null)
fi
if [[ ! -f freesasa/src/libfreesasa.a ]]; then
  git clone -q --depth 1 https://github.com/mittinatten/freesasa.git
  (cd freesasa && git submodule update --init -q && autoreconf -i >/dev/null 2>&1 \
     && ./configure --disable-json --disable-xml CFLAGS="-O3 -march=native" >/dev/null && make -j"$(nproc)" >/dev/null)
fi
gcc -O3 -march=native -fopenmp -Iparasail -Iparasail/build "$ROOT/scripts/external/pairs.c" parasail/build/libparasail.a -lm -o pairs
gcc -O3 -march=native -Ifreesasa/src "$ROOT/scripts/external/sasa.c" freesasa/src/libfreesasa.a -lm -lpthread -o sasa

SETS="random:200:100:500 fasta:data/PF00104_small.fasta random:1000:200:500 random:2000:100:400"
echo "== ours (GPU) on the alignment sets"
(cd "$ROOT" && scripts/run.sh -Dbench.skipcpu=true -Dbench.dump="$X/sets" org.biojava.tornado.bench.SwBench $SETS 2>/dev/null | grep -E "GLOBAL|LOCAL")
echo "== parasail ($(nproc) threads), same sets"
for s in $SETS; do
  tag=$(echo "$s" | sed 's/[^A-Za-z0-9]\+/_/g')
  for t in GLOBAL LOCAL; do
    m=$([[ $t == GLOBAL ]] && echo nw || echo sw)
    ./pairs "sets/$tag.fasta" $m 11 1 "sets/$tag.$t.parasail"
    cmp -s "sets/$tag.$t.scores" "sets/$tag.$t.parasail" && echo "   scores identical to ours" || echo "   SCORES DIFFER"
  done
done
echo "== ours (GPU) on the ASA structures"
(cd "$ROOT" && scripts/run.sh -Dbench.dump="$X/asa" org.biojava.tornado.bench.AsaBench 1SMT 4HHB 1AON 4V6X 3J3Q 2>/dev/null | grep -vE "cold")
echo "== FreeSASA (16 threads, its maximum for Shrake-Rupley), same atoms and radii"
for id in 1SMT 4HHB 1AON 4V6X 3J3Q; do ./sasa "asa/$id.xyzr" 16 3; done
