# 05_model.R: the Heffernan-Tawn fit shared by 05a_fit.R and 05b_simulate_chunk.R. Sourced, not run on its own,
# so both scripts fit the identical model from the same event peaks.
#
# Reads:   config.yaml (glofas.*, qc.zero_discharge_m3s, margins.*, dependence.*, simulation.*),
#          files.event_peaks
# Defines: cfg, f, sim, peaks, data (event peaks of the kept sites), zero_share, d, prior, fit (mexList),
#          m (fitted margins), margins (Laplace margins), params, limit, n_pool, n_chunks, to_p(), long()
#
# Method:
# - Stop if texmex is not 2.4.9.
# - Drop sites whose share of zero event peaks exceeds margins.max_zero_share
#   (ties at zero break the transform); print which sites were dropped.
# - mex(data, which = i, mqu = margins.threshold_quantile, dqu = dependence.threshold_quantile,
#   priorParameters) for each site i, as mexAll does but with the GPD prior (mexAll does
#   not pass it); other mex() defaults: Laplace margins, constrain = TRUE. The prior is
#   Gaussian on (phi = log sigma, xi), centred on 0, with effective sds
#   margins.phi_prior_sd and margins.xi_prior_sd. texmex 2.4.9 adds
#   (par - centre)' cov^-1 (par - centre) to -log L, without the 1/2 of a Gaussian
#   log-density, so cov = diag(2 * sd^2). No random numbers are used, so the fit is the
#   same on every machine.
# - Pool size n_pool = round(simulation.pool_factor * lambda * simulation.n_years), with
#   lambda = rows of event_peaks (historical events) / record years; it is simulated in
#   n_chunks = ceiling(n_pool / simulation.chunk_events) chunks, the last taking the remainder.
# - to_p(): per-event non-exceedance probabilities with mexTransform (the transform
#   mexMonteCarlo itself uses), then Laplace to uniform with the fit's q2p:
#   p = 0.5 * exp(y) if y < 0, else 1 - 0.5 * exp(-y). mexTransform is internal in texmex
#   2.4.9; new data go in as a copy of the fitted margins with $data replaced, against the
#   fitted margins as reference. A value above a site's GPD upper limit silently corrupts
#   that site's whole column, so stop if any value is above its limit, or if any p is not finite.

library(texmex)
library(arrow)
library(yaml)

cfg <- read_yaml("config.yaml")
f <- cfg$files
sim <- cfg$simulation
if (packageVersion("texmex") != "2.4.9") stop("texmex ", packageVersion("texmex"), " installed; 05 is written for 2.4.9")

peaks <- as.data.frame(read_parquet(f$event_peaks))
data <- peaks[, setdiff(names(peaks), "event_id")]
zero_share <- colMeans(data < cfg$qc$zero_discharge_m3s)
drop <- zero_share > cfg$margins$max_zero_share
cat(sum(drop), "sites dropped (zero share of event peaks >", cfg$margins$max_zero_share, "):", names(data)[drop], "\n")
data <- data[, !drop]
zero_share <- zero_share[!drop]
d <- ncol(data)

t0 <- Sys.time()
prior <- list(c(0, 0), diag(2 * c(cfg$margins$phi_prior_sd, cfg$margins$xi_prior_sd)^2))
prior <- setNames(rep(list(prior), d), names(data))
fit <- lapply(seq_len(d), function(i) mex(data, which = i, mqu = cfg$margins$threshold_quantile,
                                          dqu = cfg$dependence$threshold_quantile, priorParameters = prior))
names(fit) <- names(data)
oldClass(fit) <- "mexList"
cat("mex on", d, "sites and", nrow(data), "events, xi prior effective sd", cfg$margins$xi_prior_sd, ":",
    format(Sys.time() - t0, digits = 3), "\n")

m <- fit[[1]]$margins
margins <- fit[[1]]$dependence$margins
co <- coef(m)
params <- data.frame(site_id = as.integer(names(data)), threshold = co[1, ], sigma = co[3, ], xi = co[4, ],
                     zero_share = zero_share)
limit <- ifelse(params$xi < 0, params$threshold - params$sigma / params$xi, Inf)

n_pool <- round(sim$pool_factor * nrow(peaks) / (cfg$glofas$end_year - cfg$glofas$start_year + 1) * sim$n_years)
n_chunks <- ceiling(n_pool / sim$chunk_events)

ref <- m
ref$transData <- lapply(seq_len(d), function(i) m$data[, i])
to_p <- function(x, label) {
  # Non-exceedance probability of each value under the fitted margins.
  above <- sweep(as.matrix(x), 2, limit, ">")
  if (any(above)) stop(label, ": ", sum(above), " values above the GPD upper limit at sites ",
                       paste(names(data)[colSums(above) > 0], collapse = ", "))
  new <- m
  new$data <- x
  p <- margins$q2p(texmex:::mexTransform(new, margins = margins, r = ref, method = "mixture")$transformed)
  if (any(!is.finite(p))) stop(label, ": non-finite p_event")
  p
}
long <- function(p, ids) data.frame(event_id = rep(ids, times = d), site_id = rep(as.integer(names(data)), each = length(ids)),
                                    p_event = as.vector(p))
