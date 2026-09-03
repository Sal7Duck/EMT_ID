#!/usr/bin/env bash
# ssGSEA is memory-hungry; on a small machine the python process can be OOM-killed
# part way through. It checkpoints every chunk, so this just restarts it until the
# output CSV appears. Run from the repo root:  bash scores/scc_12tumour/run_ks_ssgsea.sh
cd "$(dirname "$0")/../.."
OBJ="${1:-$HOME/mnt/Sampoorn/sampoorn_figures/12_emt_scc_ID_score_HPCS_ssGSEA.h5ad}"
OUT="scores/scc_12tumour/KS_ssGSEA.csv"
for i in $(seq 1 40); do
  [ -f "$OUT" ] && break
  echo "--- attempt $i ---"
  python3 src/compute_ks_ssgsea.py "$OBJ" "$OUT" || true
done
[ -f "$OUT" ] && echo "DONE: $OUT" || echo "still not finished after 40 attempts"
