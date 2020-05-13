get_curve_properties_for_mu = function(reconstructed_model_fit, fd_type) {

# Description: For a fitted model component for "mu" that has been reconstructed on a regular grid of density ranging from zero to some positive
#              value, this function computes approximate values for some useful properties of the curve in this range.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Set the default output parameter values
curve_properties = list(q_0 = -1.0, v_ff = -1.0, dvdk_0 = -1.0, k_crit = -1.0, k_vmax = -1.0, q_cap = -1.0, v_max = -1.0, n_peaks = -1, k_jam = -1.0, v_bw = -1.0)


#######
#reconstructed_model_fit$mu[1:971] = 0.0
#reconstructed_model_fit$mu[10:21] = 1.0
#reconstructed_model_fit$mu[30:43] = 1.0
#fd_type = 'Speed.Density'
#######


# Compute some useful quantities
ngrid = nrow(reconstructed_model_fit)
grid_density_step = reconstructed_model_fit$V2[2]

# If the fitted model component for "mu" corresponds to flow
if (fd_type == 'Flow.Density') {

  # Determine the flow at zero density
  curve_properties$q_0 = reconstructed_model_fit$mu[1]

  # Estimate the free-flow speed
  curve_properties$v_ff = reconstructed_model_fit$mu[2]/grid_density_step

# If the fitted model component for "mu" corresponds to speed
} else if (fd_type == 'Speed.Density') {

  # Determine the free-flow speed
  curve_properties$v_ff = reconstructed_model_fit$mu[1]

  # Estimate the gradient of the speed (with respect to density) at zero density
  curve_properties$dvdk_0 = (reconstructed_model_fit$mu[2] - reconstructed_model_fit$mu[1])/grid_density_step
}

# Determine the grid index where the "mu" curve first becomes non-negative (if at all), and then determine the grid index at one step before the
# "mu" curve becomes negative again (if at all)
index_mu_first_nonneg = -1
index_mu_last_nonneg = -1
for (i in 1:ngrid) {
  if (reconstructed_model_fit$mu[i] < 0.0) { next }
  index_mu_first_nonneg = i
  index_mu_last_nonneg = ngrid
  for (j in index_mu_first_nonneg:ngrid) {
    if (reconstructed_model_fit$mu[j] >= 0.0) { next }
    index_mu_last_nonneg = j - 1
    break
  }
  break
}

# If all of the values in the "mu" curve are negative, then finish
if (index_mu_first_nonneg == -1) {
  return(curve_properties)
}

# In the first run of non-negative numbers in the "mu" curve, determine the grid index where the "mu" curve first becomes positive (if at all),
# and then determine the grid index at one step before the "mu" curve becomes non-positive again (if at all)
index_mu_first_pos = -1
index_mu_last_pos = -1
for (i in index_mu_first_nonneg:index_mu_last_nonneg) {
  if (reconstructed_model_fit$mu[i] <= 0.0) { next }
  index_mu_first_pos = i
  index_mu_last_pos = index_mu_last_nonneg
  for (j in index_mu_first_pos:index_mu_last_nonneg) {
    if (reconstructed_model_fit$mu[j] > 0.0) { next }
    index_mu_last_pos = j - 1
    break
  }
  break
}

# If all of the values in the first run of non-negative numbers in the "mu" curve are zero, then finish
if (index_mu_first_pos == -1) {
  return(curve_properties)
}


#### ABOVE FULLY READ AND TESTED



cat('\n')
cat(index_mu_first_nonneg, index_mu_last_nonneg, '\n')
cat(index_mu_first_pos, index_mu_last_pos, '\n')
q(save = 'no', status = 1)



cat('\n')
cat('HELLO', '\n')

cat('q_0    ', curve_properties$q_0, '\n')
cat('v_ff   ', curve_properties$v_ff, '\n')
cat('dvdk_0 ', curve_properties$dvdk_0, '\n')
cat('k_crit ', curve_properties$k_crit, '\n')
cat('k_vmax ', curve_properties$k_vmax, '\n')
cat('q_cap  ', curve_properties$q_cap, '\n')
cat('v_max  ', curve_properties$v_max, '\n')
cat('n_peaks', curve_properties$n_peaks, '\n')
cat('k_jam  ', curve_properties$k_jam, '\n')
cat('v_bw   ', curve_properties$v_bw, '\n')

q(save = 'no', status = 1)


  # Estimate the flow capacity
#  k_crit_sub = which.max(reconstructed_model_fit$mu)
#  k_crit = reconstructed_model_fit$V2[k_crit_sub]
#  q_cap = reconstructed_model_fit$mu[k_crit_sub]

  #
#  pp = findpeaks(reconstructed_model_fit$mu)



a = 1

# Return the ....
return(a)
}
