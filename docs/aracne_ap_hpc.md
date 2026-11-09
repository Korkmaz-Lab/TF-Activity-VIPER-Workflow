# Running ARACNe-AP (notes and SLURM templates)

Short practical notes for building a network with [ARACNe-AP](https://github.com/califano-lab/ARACNe-AP) (Lachmann et al. 2016). Option names and behaviour are the tool's own; check them against the version you use (`java -jar Aracne.jar` prints usage). The scripts below are **templates**: replace everything in `<angle brackets>` and adapt to your cluster.

## Inputs

| File | Format |
|---|---|
| `network_input.tsv` | Gene × sample, tab-delimited. First column = gene symbol, first row = sample names. MI is rank-based, so VST, log2(TPM+1) or log2(FPKM+0.1) all work; state which one you used. One cohort/context per network (typically ≥ 100 samples). |
| `regulators.txt` | One gene symbol per line (e.g. a curated TF list). Every regulator must be a row of `network_input.tsv`. |

Do not open or save gene tables in Excel (SEPT/MARCH symbols turn into dates). Keep `network_input.tsv` unchanged after the run: the same matrix is needed later for `aracne2regulon()` (see the main README).

## The three steps

| Step | What it does | Job type |
|---|---|---|
| 1. MI threshold | `--calculateThreshold` writes `miThreshold_p<pvalue>_samples<N>.txt`. It depends on the number of samples, the p-value and the seed, not on the expression values or the regulator list, so it is computed once per sample number and reused. | one job |
| 2. Bootstraps | One network per `--seed`. The threshold file must be in the output directory. | array, one task per seed |
| 3. Consolidate | `--consolidate` merges every `bootstrapNetwork_*` file in a directory into `network.txt`. | one job after the array |

Settings: MI p-value `1E-8`, DPI tolerance `0` and Bonferroni consolidation are ARACNe-AP defaults. **100 bootstraps is common practice, not a default.** Write the values you used next to the results.

`network.txt` has a header and four columns: `Regulator`, `Target`, `MI`, `pvalue`.

## Java 8 and Singularity

On the lab's cluster `Aracne.jar` did not run directly with the Java modules installed there, so it is run inside a Java 8 container image (`java8.sif`; any image with a Java 8 runtime works):

```bash
singularity exec <java8_sandbox> java -Xmx<N>G -jar Aracne.jar ...
```

Two cluster-specific points worth checking on yours:

- If FUSE mounts are not allowed on compute nodes, unpack the image once into a sandbox directory (`singularity build --sandbox <sandbox_dir> java8.sif`) and run every task from that directory instead of from the `.sif`.
- `SINGULARITY_CACHEDIR` must point to a writable directory (e.g. `$HOME/.singularity/cache`).

Record the `md5sum` of `Aracne.jar`, the image and the input files in the job logs (the templates do this). Also note where you obtained the jar (release or source commit).

## Resources (lab experience, will vary with your data)

- Threshold job: minutes, ~12 GiB for ~3,000 samples.
- Bootstrap task: ~1.5 k regulators × ~3 k samples took 14–20 min and ~15 GiB; ~7 k regulators took 1.7–3.2 h and ≤ 36 GiB. MI time grows with the number of regulators, DPI time roughly with its square. Set the Java heap (`-Xmx`) below the SLURM `--mem`.
- Run one pilot bootstrap first to measure time and memory. Do not mix pilot bootstraps into the final 100.
- If one array task fails, the consolidation job (`afterok`) never starts: rerun the failed seeds, then submit consolidation again.

## Templates

Run from a project directory; paths are relative to it. Layout used below:

```
software/Aracne.jar   software/java8_sandbox/   input/network_input.tsv   input/regulators.txt
network/threshold/    network/bootstraps/       network/network.txt
```

**`setup_sandbox.sh`** (once)

```bash
#!/bin/bash
#SBATCH --job-name=aracne_setup
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --partition=<partition>
set -euo pipefail
module load singularity                      # site-specific
export SINGULARITY_CACHEDIR=$HOME/.singularity/cache
mkdir -p "$SINGULARITY_CACHEDIR"
singularity build --sandbox software/java8_sandbox software/java8.sif
singularity exec software/java8_sandbox java -version
```

**`threshold.sh`** (step 1)

```bash
#!/bin/bash
#SBATCH --job-name=aracne_thr
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=<partition>
set -euo pipefail
module load singularity
export SINGULARITY_CACHEDIR=$HOME/.singularity/cache
SW=software; EXPR=input/network_input.tsv; REG=input/regulators.txt; OUT=network/threshold
mkdir -p "$OUT"
md5sum $SW/Aracne.jar $EXPR $REG
singularity exec $SW/java8_sandbox java -Xmx24G -jar $SW/Aracne.jar \
  -e $EXPR -o "$OUT" --tfs $REG --pvalue 1E-8 --seed 1 --calculateThreshold
ls "$OUT"/miThreshold_p1E-8_samples*.txt
```

**`bootstrap_array.sh`** (step 2; submit with `--array=1-100`, the task ID is the seed)

```bash
#!/bin/bash
#SBATCH --job-name=aracne_boot
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=<hh:mm:ss>
#SBATCH --partition=<partition>
#SBATCH --array=1-100
set -euo pipefail
module load singularity
export SINGULARITY_CACHEDIR=$HOME/.singularity/cache
SW=software; EXPR=input/network_input.tsv; REG=input/regulators.txt
SEED=$SLURM_ARRAY_TASK_ID; TAG=$(printf 'seed%03d' "$SEED")
OUT=network/bootstraps; mkdir -p "$OUT"
WORK=$(mktemp -d "${TMPDIR:-/tmp}/aracne_${SLURM_JOB_ID}_XXXX")   # node-local scratch
trap 'rm -rf "$WORK"' EXIT
cp network/threshold/miThreshold_p1E-8_samples*.txt "$WORK"/      # the jar looks for it in -o
singularity exec $SW/java8_sandbox java -Xmx32G -jar $SW/Aracne.jar \
  -e $EXPR -o "$WORK" --tfs $REG --pvalue 1E-8 --seed "$SEED" --threads "$SLURM_CPUS_PER_TASK"
[ "$(ls "$WORK"/bootstrapNetwork_*.txt | wc -l)" -eq 1 ] || { echo "expected 1 bootstrap file"; exit 1; }
cp "$WORK"/bootstrapNetwork_*.txt "$OUT/bootstrapNetwork_${TAG}.txt"
```

**`consolidate.sh`** (step 3)

```bash
#!/bin/bash
#SBATCH --job-name=aracne_cons
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=<partition>
set -euo pipefail
module load singularity
export SINGULARITY_CACHEDIR=$HOME/.singularity/cache
SW=software; D=network
[ "$(ls $D/bootstraps/bootstrapNetwork_seed*.txt | wc -l)" -eq 100 ] || { echo "expected 100 bootstraps"; exit 1; }
md5sum $D/bootstraps/bootstrapNetwork_seed*.txt > $D/bootstrap_md5.txt
singularity exec $SW/java8_sandbox java -Xmx24G -jar $SW/Aracne.jar -o $D/bootstraps --consolidate
mv $D/bootstraps/network.txt $D/network.txt
echo "edges: $(($(wc -l < $D/network.txt) - 1))  regulators: $(tail -n +2 $D/network.txt | cut -f1 | sort -u | wc -l)"
```

**Submitting** (login node)

```bash
T=$(sbatch --parsable threshold.sh)
A=$(sbatch --parsable --dependency=afterok:$T bootstrap_array.sh)
sbatch --dependency=afterok:$A consolidate.sh
# rerun failed seeds only, e.g.:  sbatch --array=37,52 bootstrap_array.sh   (then resubmit consolidate.sh)
```

Next: load `network/network.txt` in R as described in the main README.
