# 05a_fit.R: fit the Heffernan-Tawn model with texmex (on the PC); the simulation runs on the cluster (05b).
#
# Reads:  config.yaml, files.event_peaks (via 05_model.R), files.sites
# Writes: files.margins_params (site_id, threshold, sigma, xi, zero_share),
#         files.dependence_params (cond_site, site_id, a, b): Heffernan-Tawn alpha and beta,
#         files.hist_probs (event_id, site_id, p_event) for the historical events,
#         diagnostic PNGs for the required (hotspot) sites in files.figures
#
# Method:
# - Fit as in 05_model.R. Print xi and, where xi < 0, the GPD upper limit at each site.
# - hist_probs: the historical event peaks through to_p().
# - Diagnostics, conditioning on each required site, at a fixed size: gpdRangeFit
#   (thresholds between the margins.range_quantiles quantiles; texmex's default 0.05-0.95
#   range fails at Budalangi), and plot() and mexRangeFit of a fit to the required sites'
#   columns only. Heffernan-Tawn dependence is fitted pair by pair, so that fit's alpha and
#   beta equal the full fit's; print the largest difference, and stop above 1e-6.
# - Print n_pool, the number of chunks n_chunks, the sbatch command for the cluster and the chunk folder.

source("scripts/05_model.R")

dep <- do.call(rbind, lapply(fit, function(x) {
  co <- coef(x$dependence)
  data.frame(cond_site = as.integer(names(x$margins$data)[x$dependence$which]), site_id = as.integer(colnames(co)),
             a = co["a", ], b = co["b", ])
}))
dir.create(dirname(f$dependence_params), showWarnings = FALSE, recursive = TRUE)
write_parquet(dep, f$dependence_params)
write_parquet(params, f$margins_params)
cat("\nxi by site, with the GPD upper limit where xi < 0 (m3/s):\n")
tab <- data.frame(site_id = params$site_id, xi = round(params$xi, 3), threshold = round(params$threshold, 1),
                  upper_limit = round(ifelse(is.finite(limit), limit, NA), 1), record_max = round(sapply(data, max), 1))
print(tab[order(tab$xi), ], row.names = FALSE)
cat(sum(params$xi < 0), "of", d, "sites with xi < 0\n")

hist <- long(to_p(data, "historical"), peaks$event_id)
write_parquet(hist, f$hist_probs)
cat("wrote", f$hist_probs, "(", nrow(hist), "rows )\n")

# Diagnostics, conditioning on each required site, on the required sites' columns.
sites <- as.data.frame(read_parquet(f$sites))
required <- sites[sites$name != "" & sites$site_id %in% params$site_id, ]
cols <- as.character(required$site_id)
sub <- data[, cols]
dir.create(f$figures, showWarnings = FALSE, recursive = TRUE)
worst <- 0
for (k in seq_along(cols)) {
  stem <- file.path(f$figures, paste0("mex_", required$name[k]))
  sfit <- mex(sub, which = k, mqu = cfg$margins$threshold_quantile, dqu = cfg$dependence$threshold_quantile,
              priorParameters = prior[cols])
  full <- coef(fit[[cols[k]]]$dependence)[, setdiff(cols, cols[k]), drop = FALSE]
  worst <- max(worst, abs(coef(sfit$dependence)[c("a", "b"), colnames(full)] - full[c("a", "b"), ]))
  png(paste0(stem, "_fit.png"), width = 1200, height = 320 * (length(cols) - 1))
  par(mfrow = c(length(cols) - 1, 3))
  plot(sfit)
  dev.off()
  png(paste0(stem, "_gpd_range.png"), width = 900, height = 450)
  par(mfrow = c(1, 2))
  q <- quantile(data[, cols[k]], cfg$margins$range_quantiles)
  plot(gpdRangeFit(data[, cols[k]], umin = q[1], umax = q[2]))
  dev.off()
  png(paste0(stem, "_mex_range.png"), width = 900, height = 320 * (length(cols) - 1))
  par(mfrow = c(length(cols) - 1, 2))
  plot(mexRangeFit(sfit$margins, which = k, trace = Inf))
  dev.off()
  cat("diagnostics for", required$name[k], "written to", paste0(stem, "_*.png"), "\n")
}
cat("largest |difference| in alpha, beta between the required-site fit and the full fit:", signif(worst, 3), "\n")
if (worst > 1e-6) stop("required-site fit differs from the full fit")

cat("\npool:", n_pool, "events (pool_factor", sim$pool_factor, "x lambda x", sim$n_years, "years) in", n_chunks,
    "chunks of", sim$chunk_events, "\n")
cat("on the cluster: sbatch --array=1-", n_chunks, "%100 cluster/run_chunks.sbatch\n", sep = "")
cat("chunks go to", chunk_dir, "\n")
