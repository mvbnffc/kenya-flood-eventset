# 05_simulate.R: fit the Heffernan-Tawn model with texmex and simulate events.
#
# Reads:  config.yaml (glofas.*, qc.zero_discharge_m3s, margins.*, dependence.*, simulation.*),
#         files.event_peaks, files.events, files.sites
# Writes: files.mex_fit (fitted mexList, .rds),
#         files.margins_params (site_id, threshold, sigma, xi, zero_share),
#         files.hist_probs (event_id, site_id, p_event) for the historical events,
#         files.sim_probs (event_id, sim_year, source, site_id, p_event),
#         diagnostic PNGs in files.figures
#
# Method:
# - Drop sites whose share of zero event peaks exceeds margins.max_zero_share
#   (ties at zero break the transform); print which sites were dropped.
# - mexAll(data, mqu = margins.threshold_quantile, dqu = rep(dependence.threshold_quantile, d)).
#   mexAll only takes data, mqu and dqu, so it uses mex() defaults: Laplace margins
#   and constrain = TRUE (checked in texmex 2.4.9). Print xi and, where xi < 0, the
#   GPD upper limit at each site.
# - Diagnostics for the required sites as conditioning sites: plot() of the mex fit,
#   gpdRangeFit (thresholds between the margins.range_quantiles quantiles; texmex's
#   default 0.05-0.95 range fails at Budalangi) and mexRangeFit, saved as PNGs.
# - Annual event counts (simulation.count_model = empirical): count the historical
#   events in each record year (zeros included), then draw simulation.n_years counts
#   from those with replacement, as in Quinn et al. (2019). This keeps the clustering
#   of events in wet years without fitting a count model. Stop if count_model is
#   anything else.
# - mexMonteCarlo(sum(counts), fit, mult = simulation.mult). Its MCsample is on
#   the original discharge scale. It asks predict() for sum(counts) * d * mult
#   samples per conditioning site, so run time scales with mult; simulation.mult is
#   the lowest value with no NA rows in a 500-event test, plus one. Stop on NA rows
#   (too few kept samples). For the 10,000-year run, parallel chunks of years with
#   parallel::parLapply and clusterSetRNGStream are planned (not written yet).
# - Convert simulated and historical peaks to per-event non-exceedance probabilities
#   p_event with mexTransform (the transform mexMonteCarlo itself uses), then Laplace
#   to uniform with the fit's q2p: p = 0.5 * exp(y) if y < 0, else 1 - 0.5 * exp(-y).
#   mexTransform is internal in texmex 2.4.9; new data go in as a copy of the fitted
#   margins with $data replaced, against the fitted margins as reference. A value
#   above a site's GPD upper limit silently corrupts that site's whole column, so
#   stop if any value is above its limit. Stop if texmex is not 2.4.9.
# - source = "model" where whichMaxAboveThresh is TRUE, else "resampled".
# - Assign simulated events to years in order using the counts.

library(texmex)
library(arrow)
library(yaml)

cfg <- read_yaml("config.yaml")
set.seed(cfg$simulation$seed)
f <- cfg$files
sim <- cfg$simulation
if (packageVersion("texmex") != "2.4.9") stop("texmex ", packageVersion("texmex"), " installed; 05 is written for 2.4.9")
if (sim$count_model != "empirical") stop("simulation.count_model must be empirical, not ", sim$count_model)

peaks <- as.data.frame(read_parquet(f$event_peaks))
events <- as.data.frame(read_parquet(f$events))
if (!identical(peaks$event_id, events$event_id)) stop("event_peaks and events have different event_ids")
data <- peaks[, setdiff(names(peaks), "event_id")]
zero_share <- colMeans(data < cfg$qc$zero_discharge_m3s)
drop <- zero_share > cfg$margins$max_zero_share
cat(sum(drop), "sites dropped (zero share of event peaks >", cfg$margins$max_zero_share, "):",
    names(data)[drop], "\n")
data <- data[, !drop]
zero_share <- zero_share[!drop]
d <- ncol(data)

t0 <- Sys.time()
fit <- mexAll(data, mqu = cfg$margins$threshold_quantile, dqu = rep(cfg$dependence$threshold_quantile, d))
cat("mexAll on", d, "sites and", nrow(data), "events:", format(Sys.time() - t0, digits = 3), "\n")
dir.create(dirname(f$mex_fit), showWarnings = FALSE, recursive = TRUE)
saveRDS(fit, f$mex_fit)

m <- fit[[1]]$margins
margins <- fit[[1]]$dependence$margins
co <- coef(m)
params <- data.frame(site_id = as.integer(names(data)), threshold = co[1, ], sigma = co[3, ], xi = co[4, ],
                     zero_share = zero_share)
write_parquet(params, f$margins_params)
limit <- ifelse(params$xi < 0, params$threshold - params$sigma / params$xi, Inf)
cat("\nxi by site, with the GPD upper limit where xi < 0 (m3/s):\n")
tab <- data.frame(site_id = params$site_id, xi = round(params$xi, 3), threshold = round(params$threshold, 1),
                  upper_limit = round(ifelse(is.finite(limit), limit, NA), 1), record_max = round(sapply(data, max), 1))
print(tab[order(tab$xi), ], row.names = FALSE)

# Diagnostics, conditioning on each required site.
sites <- as.data.frame(read_parquet(f$sites))
required <- sites[sites$name != "" & sites$site_id %in% params$site_id, ]
dir.create(f$figures, showWarnings = FALSE, recursive = TRUE)
for (k in seq_len(nrow(required))) {
  i <- match(as.character(required$site_id[k]), names(data))
  stem <- file.path(f$figures, paste0("mex_", required$name[k]))
  png(paste0(stem, "_fit.png"), width = 1200, height = 320 * (d - 1))
  par(mfrow = c(d - 1, 3))
  plot(fit[[i]])
  dev.off()
  png(paste0(stem, "_gpd_range.png"), width = 900, height = 450)
  par(mfrow = c(1, 2))
  q <- quantile(data[, i], cfg$margins$range_quantiles)
  plot(gpdRangeFit(data[, i], umin = q[1], umax = q[2]))
  dev.off()
  png(paste0(stem, "_mex_range.png"), width = 900, height = 320 * (d - 1))
  par(mfrow = c(d - 1, 2))
  plot(mexRangeFit(m, which = i, trace = Inf))
  dev.off()
  cat("diagnostics for", required$name[k], "written to", paste0(stem, "_*.png"), "\n")
}

# Annual counts, resampled from the historical record years.
years <- cfg$glofas$start_year:cfg$glofas$end_year
annual <- as.vector(table(factor(as.integer(format(events$start, "%Y", tz = "UTC")), levels = years)))
counts <- sample(annual, sim$n_years, replace = TRUE)
cat("\nlambda =", nrow(events), "events /", length(years), "years =", round(nrow(events) / length(years), 2),
    "; historical annual counts", min(annual), "to", max(annual), "; simulated events", sum(counts), "\n")

t0 <- Sys.time()
mc <- mexMonteCarlo(sum(counts), fit, mult = sim$mult)
cat("mexMonteCarlo:", format(Sys.time() - t0, digits = 3), "\n")
if (any(!complete.cases(mc$MCsample))) stop(sum(!complete.cases(mc$MCsample)), " NA rows in MCsample; raise simulation.mult")

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

hist <- long(to_p(data, "historical"), events$event_id)
write_parquet(hist, f$hist_probs)
n_sim <- sum(counts)
p_sim <- long(to_p(mc$MCsample, "simulated"), seq_len(n_sim))
p_sim$sim_year <- rep(rep(seq_len(sim$n_years), counts), times = d)
p_sim$source <- rep(ifelse(mc$whichMaxAboveThresh, "model", "resampled"), times = d)
write_parquet(p_sim[, c("event_id", "sim_year", "source", "site_id", "p_event")], f$sim_probs)
cat("wrote", f$hist_probs, "(", nrow(hist), "rows ) and", f$sim_probs, "(", nrow(p_sim), "rows ):",
    n_sim, "simulated events in", sim$n_years, "years,", sum(mc$whichMaxAboveThresh), "from the model\n")
