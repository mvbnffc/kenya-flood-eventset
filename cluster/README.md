# Running the simulation (05b) on the SLURM cluster

05a (fit), 05c (combine), 06 and 07 run on the PC. Only 05b runs on the cluster: one array task per chunk of `simulation.chunk_events` events. Each task refits the model from `event_peaks.parquet` and uses its own random-number stream, so any chunk can be rerun alone and gives the same output.

The cluster needs system R (`/usr/bin/Rscript`) with texmex 2.4.9, arrow, mvtnorm and yaml in your user library, and a clone of this repo. Below, `CLUSTER` stands for `user@cluster-host` and `REPO` for the path of the clone on the cluster.

## On the PC, before the cluster run

1. Run 02-04, then `Rscript scripts/05a_fit.R`. Its last line is the sbatch command with the chunk count N, e.g. `sbatch --array=1-38%100 cluster/run_chunks.sbatch`.
2. Commit and push the config and scripts the run should use (the cluster reads `config.yaml` from git).

## On the cluster

```
cd REPO
git pull
mkdir -p data/interim      # first time only: data/ is not in git
```

From the PC (Git Bash, which has scp but not rsync), copy the event peaks over:

```
scp data/interim/event_peaks.parquet CLUSTER:REPO/data/interim/
```

Then, on the cluster, from `REPO`:

```
sbatch --array=1-N%100 cluster/run_chunks.sbatch
```

with N from step 1. Logs go to `cluster/logs/chunk_<chunk>_<job>.out` and `.err`; each `.out` ends with the chunk's events, time per event and peak memory. A chunk that fails writes no file; rerun it alone with `sbatch --array=<chunk> cluster/run_chunks.sbatch`.

## Back on the PC

Copy the chunks back (into `results/chunks/`) and combine them:

```
mkdir -p results
scp -r CLUSTER:REPO/results/chunks results/
python scripts/05c_combine.py
python scripts/06_event_rp.py
python scripts/07_checks.py
```

05c stops if any chunk is missing or has NA values; rerun those chunks on the cluster and copy again.
