# 05_simulate.R: fit the Heffernan-Tawn model with texmex and simulate events.
#
# Reads:  config.yaml (glofas.*, margins.*, dependence.*, simulation.*), files.event_peaks, files.events
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
#   and constrain = TRUE. Check those defaults in the installed texmex version.
#   If mexAll is too slow at full scale, copy its five lines and swap lapply for
#   parallel::mclapply.
# - Diagnostics for five representative conditioning sites: plot() of the mex fit,
#   gpdRangeFit and mexRangeFit, saved as PNGs.
# - Annual event counts (simulation.count_model = empirical): count the historical
#   events in each record year (zeros included), then draw simulation.n_years counts
#   from those with replacement, as in Quinn et al. (2019). This keeps the clustering
#   of events in wet years without fitting a count model. Stop if count_model is
#   anything else.
# - mexMonteCarlo(sum(counts), fit, mult = simulation.mult). Its MCsample is on
#   the original discharge scale.
# - Convert simulated and historical peaks to per-event non-exceedance probabilities
#   p_event with mexTransform (the transform mexMonteCarlo itself uses), then Laplace
#   to uniform: p = 0.5 * exp(y) if y < 0, else 1 - 0.5 * exp(-y). First check that
#   mexTransform accepts new data; if it does not, stop and ask.
# - source = "model" where whichMaxAboveThresh is TRUE, else "resampled".
# - Assign simulated events to years in order using the counts.
# - If mexMonteCarlo runs out of memory, simulate in chunks of years.

library(texmex)
library(arrow)
library(yaml)

cfg <- read_yaml("config.yaml")
set.seed(cfg$simulation$seed)

# To be written in an agent session.
