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

# Determine the grid index where the "mu" curve first becomes positive (if at all), and then determine the grid index at one step before the "mu"
# curve becomes non-positive again (if at all)
index_mu_first_pos = -1
index_mu_last_pos = -1
for (i in 1:ngrid) {
  if (reconstructed_model_fit$mu[i] <= 0.0) { next }
  index_mu_first_pos = i
  index_mu_last_pos = ngrid
  for (j in index_mu_first_pos:ngrid) {
    if (reconstructed_model_fit$mu[j] > 0.0) { next }
    index_mu_last_pos = j - 1
    break
  }
  break
}

# If all of the values in the "mu" curve are non-positive, then finish
if (index_mu_first_pos == -1) {
  return(curve_properties)
}

# Create a temporary vector that starts with a single zero value, followed by the first run of positive numbers in the "mu" curve, and ends with
# another single zero value
npos_p2 = index_mu_last_pos - index_mu_first_pos + 3
tmpvec = double(length = npos_p2)
tmpvec[2:(npos_p2 - 1)] = reconstructed_model_fit$mu[index_mu_first_pos:index_mu_last_pos]

# Find the number of peaks in the first run of positive numbers in the "mu" curve while also finding the highest peak and its corresponding density.
# In the case that two or more peaks are all equally the highest peak, then the peak with the smallest density is selected as the highest peak.
curr_i = 2
n_peaks = 0
max_peak_value = -1.0
max_peak_density = -1.0
while (curr_i < npos_p2) {

  # Extract the previous, current, and next values from the temporary vector
  prev_tmpvec = tmpvec[curr_i - 1]
  curr_tmpvec = tmpvec[curr_i]
  next_tmpvec = tmpvec[curr_i + 1]

  # If the previous value in the temporary vector is greater than the current value, then there is no peak at the current value
  if (prev_tmpvec > curr_tmpvec) {

    # If the next value in the temporary vector is not equal to the current value, then move on to the next value
    if (next_tmpvec != curr_tmpvec) {
      curr_i = curr_i + 1
      next
    }

    # If the next value in the temporary vector is equal to the current value, then move on to the next value that is not equal to the current
    # value
    while (tmpvec[curr_i] == curr_tmpvec) { curr_i = curr_i + 1 }
    next

  # If the previous value in the temporary vector is less than the current value, then it is possible that there is a peak at, or near, the current
  # value (N.B: Due to the way that this algorithm is designed, it is impossible that the previous value in the temporary vector is equal to the
  # current value)
  } else {

    # If the next value in the temporary vector is greater than the current value, then there is no peak at the current value. Move on to the next
    # value.
    if (next_tmpvec > curr_tmpvec) {
      curr_i = curr_i + 1
      next
    }

    # If the next value in the temporary vector is less than the current value, then there is a peak at the current value. Count the peak and update
    # the highest peak if necessary. Move on to the next value.
    if (next_tmpvec < curr_tmpvec) {
      n_peaks = n_peaks + 1
      if (curr_tmpvec > max_peak_value) {
        max_peak_value = curr_tmpvec
        max_peak_density = reconstructed_model_fit$V2[curr_i + index_mu_first_pos - 2]
      }
      curr_i = curr_i + 1
      next
    }

    # If the next value in the temporary vector is equal to the current value, then find the next value that is not equal to the current value
    start_i = curr_i
    while (tmpvec[curr_i] == curr_tmpvec) { curr_i = curr_i + 1 }

    # If the next value that is not equal to the current value is greater than the current value, then there is no peak near to the current value.
    # Move on to the next value that is not equal to the current value.
    if (tmpvec[curr_i] > curr_tmpvec) { next }

    # If the next value that is not equal to the current value is less than the current value, then there is a peak near to the current value.
    # Count the peak and update the highest peak if necessary. Move on to the next value that is not equal to the current value.
    n_peaks = n_peaks + 1
    if (curr_tmpvec > max_peak_value) {
      max_peak_value = curr_tmpvec
      max_peak_density = 0.5*(reconstructed_model_fit$V2[start_i + index_mu_first_pos - 2] + reconstructed_model_fit$V2[curr_i + index_mu_first_pos - 3])
    }
  }
}


#### ABOVE FULLY READ AND TESTED


cat('\n')
cat(index_mu_first_pos, index_mu_last_pos, '\n')

cat('\n')
cat(reconstructed_model_fit$V2[index_mu_first_pos:index_mu_last_pos], '\n')
cat('\n')
cat(tmpvec, '\n')


cat('\n')
cat(n_peaks, '\n')
cat(max_peak_value, '\n')
cat(max_peak_density, '\n')

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
