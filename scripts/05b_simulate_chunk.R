# 05b_simulate_chunk.R: simulate one chunk of the event pool (on the cluster, one SLURM array task per chunk).
#
# Usage:  Rscript scripts/05b_simulate_chunk.R <chunk_id>     (1 to n_chunks; 05a prints n_chunks)
# Reads:  config.yaml, files.event_peaks (via 05_model.R, which refits the model)
# Writes: <chunk_dir>/chunk_NNNN.parquet (event_id, source, site_id, p_event); chunk_dir is files.chunks
#         with {seed} and {n_years} filled from simulation.*, so runs with other settings can't mix
#
# Method:
# - Refit as in 05_model.R (no random numbers, so identical to 05a's fit).
# - Random numbers: RNGkind("L'Ecuyer-CMRG"), set.seed(simulation.seed), then
#   parallel::nextRNGStream applied chunk_id times. Each chunk has its own stream, so a
#   chunk gives the same output whenever and wherever it is run.
# - mexMonteCarlo(n, fit, mult = simulation.mult) with n = simulation.chunk_events, or the
#   remainder of n_pool for the last chunk. It asks predict() for n * d * mult samples per
#   conditioning site, so run time scales with mult. Stop on NA rows (too few kept samples),
#   so a failed chunk writes no file.
# - p_event with to_p(); event_id = (chunk_id - 1) * simulation.chunk_events + 1 ... n.
#   source = "model" where whichMaxAboveThresh is TRUE, else "resampled".
# - Print the time per event and R's peak memory (gc max used).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript scripts/05b_simulate_chunk.R <chunk_id>")
id <- as.integer(args[1])

source("scripts/05_model.R")
if (is.na(id) || id < 1 || id > n_chunks) stop("chunk_id must be 1 to ", n_chunks, ", not ", args[1])
n <- if (id < n_chunks) sim$chunk_events else n_pool - (n_chunks - 1) * sim$chunk_events

RNGkind("L'Ecuyer-CMRG")
set.seed(sim$seed)
stream <- .Random.seed
for (k in seq_len(id)) stream <- parallel::nextRNGStream(stream)
assign(".Random.seed", stream, envir = globalenv())

gc(reset = TRUE)
t0 <- Sys.time()
mc <- mexMonteCarlo(n, fit, mult = sim$mult)
secs <- as.numeric(Sys.time() - t0, units = "secs")
if (any(!complete.cases(mc$MCsample))) stop(sum(!complete.cases(mc$MCsample)), " NA rows in MCsample; raise simulation.mult")

ids <- as.integer((id - 1) * sim$chunk_events + seq_len(n))
out <- long(to_p(mc$MCsample, paste("chunk", id)), ids)
out$source <- rep(ifelse(mc$whichMaxAboveThresh, "model", "resampled"), times = d)
dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
path <- file.path(chunk_dir, sprintf("chunk_%04d.parquet", id))
write_parquet(out[, c("event_id", "source", "site_id", "p_event")], path)
g <- gc()
cat(sprintf("chunk %d of %d: events %d-%d, %d from the model; mexMonteCarlo %.0f s (%.2f s per event, mult %d); R peak %.0f Mb; wrote %s\n",
            id, n_chunks, min(ids), max(ids), sum(mc$whichMaxAboveThresh), secs, secs / n, sim$mult,
            sum(g[, ncol(g)]), path))
