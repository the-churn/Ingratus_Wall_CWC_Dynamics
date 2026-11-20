# ==============================================================================
# age_model.R
# The age model (ESM Methods S3), sourced by the scripts that compute ages.
# A colony's modelled age is the time the fitted growth curve needs to take
# it from the detection floor A0 to its observed planar area A:
#
#   t(A) = integral from A0 to A of  da / (a * RGR(a)),
#   RGR(a) = asym + (r0 - asym) * exp(-exp(lrc) * a)
# ==============================================================================

# Species-level parameters (asym, r0, lrc) for each row of a matrix of
# parameter vectors: the 5,000 draws, or the fitted coefficients (a named
# vector becomes one row). M. oculata is the model's reference level.
species_params <- function(draws, species) {
  draws <- rbind(draws)
  g <- function(par) {
    draws[, paste0(par, ".(Intercept)")] +
      if (species == SPECIES[1]) 0 else draws[, paste0(par, ".fSpecies", species)]
  }
  cbind(asym = g("Asym"), r0 = g("R0"), lrc = g("lrc"))
}

# Modelled age (yr) of a colony of area `size` (cm2), for one parameter set
# p (a row of species_params()). Zero at or below A0; NA if the curve
# reaches RGR <= 0 below that size (the colony could not have grown there).
age_at_size <- function(size, p, A0) {
  if (size <= A0) return(0)
  f_integrand <- function(s) {
    rgr <- p[["asym"]] + (p[["r0"]] - p[["asym"]]) * exp(-exp(p[["lrc"]]) * s)
    if (any(rgr <= 0, na.rm = TRUE)) return(rep(NA_real_, length(s)))
    1 / (s * rgr)
  }
  tryCatch({
    val <- integrate(f_integrand, lower = A0, upper = size)$value
    if (is.nan(val) || val < 0) NA_real_ else val
  }, error = function(e) NA_real_)
}
