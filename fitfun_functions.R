################################################################################################################################################
get_npeaks_ntroughs = function(xvec, yvec, edge_value) {

# Description: For a curve in the xy-plane that is represented by a discrete set of N points with (x,y) pairs "xvec" and "yvec", this function
#              determines the number of peaks and troughs, and the coordinates of the highest peak and the lowest trough. A subtle but important
#              requirement for this algorithm to work properly is that "edge_value" must not be not equal to the first or last elements of "yvec".
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create a temporary vector that starts with a single value "edge_value", followed by the vector of numbers "yvec", and ends with another single
# value "edge_value"
N_p2 = length(yvec) + 2
tmpvec = double(length = N_p2)
tmpvec[1] = edge_value
tmpvec[2:(N_p2 - 1)] = yvec
tmpvec[N_p2] = edge_value

# Find the number of peaks and troughs in the temporary vector while also finding the highest peak and the lowest trough and their corresponding
# x-values. In the case that two or more peaks are all equally the highest peak, then the peak with the lowest x-value is selected as the highest
# peak. In the case that two or more troughs are all equally the lowest trough, then the trough with the lowest x-value is selected as the lowest
# trough.
curr_i = 2
n_peaks = 0
highest_peak_x = 0.0
highest_peak_y = 0.0
n_troughs = 0
lowest_trough_x = 0.0
lowest_trough_y = 0.0
while (curr_i < N_p2) {

  # Extract the previous, current, and next values from the temporary vector
  prev_tmpvec = tmpvec[curr_i - 1]
  curr_tmpvec = tmpvec[curr_i]
  next_tmpvec = tmpvec[curr_i + 1]

  # If the previous value in the temporary vector is greater than the current value, then there is no peak at the current value. However, it is
  # possible that there is a trough at, or near, the current value.
  if (prev_tmpvec > curr_tmpvec) {

    # If the next value in the temporary vector is greater than the current value, then there is a trough at the current value. Count the trough
    # and update the lowest trough if necessary. Move on to the next value.
    if (next_tmpvec > curr_tmpvec) {
      n_troughs = n_troughs + 1
      if (n_troughs == 1) {
        lowest_trough_x = xvec[curr_i - 1]
        lowest_trough_y = curr_tmpvec
      } else {
        if (curr_tmpvec < lowest_trough_y) {
          lowest_trough_x = xvec[curr_i - 1]
          lowest_trough_y = curr_tmpvec
        }
      }
      curr_i = curr_i + 1
      next
    }

    # If the next value in the temporary vector is less than the current value, then there is no trough at the current value. Move on to the next
    # value.
    if (next_tmpvec < curr_tmpvec) {
      curr_i = curr_i + 1
      next
    }

    # If the next value in the temporary vector is equal to the current value, then find the next value that is not equal to the current value
    start_i = curr_i
    while (tmpvec[curr_i] == curr_tmpvec) { curr_i = curr_i + 1 }

    # If the next value that is not equal to the current value is less than the current value, then there is no trough near to the current value.
    # Move on to the next value that is less than the current value.
    if (tmpvec[curr_i] < curr_tmpvec) { next }

    # If the next value that is not equal to the current value is greater than the current value, then there is a trough near to the current value.
    # Count the trough and update the lowest trough if necessary. Move on to the next value that is greater than the current value.
    n_troughs = n_troughs + 1
    if (n_troughs == 1) {
      lowest_trough_x = 0.5*(xvec[start_i - 1] + xvec[curr_i - 2])
      lowest_trough_y = curr_tmpvec
    } else {
      if (curr_tmpvec < lowest_trough_y) {
        lowest_trough_x = 0.5*(xvec[start_i - 1] + xvec[curr_i - 2])
        lowest_trough_y = curr_tmpvec
      }
    }
    next

  # If the previous value in the temporary vector is less than the current value, then there is no trough at the current value. However, it is
  # possible that there is a peak at, or near, the current value (N.B: Due to the way that this algorithm is designed, it is impossible that the
  # previous value in the temporary vector is equal to the current value)
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
      if (n_peaks == 1) {
        highest_peak_x = xvec[curr_i - 1]
        highest_peak_y = curr_tmpvec
      } else {
        if (curr_tmpvec > highest_peak_y) {
          highest_peak_x = xvec[curr_i - 1]
          highest_peak_y = curr_tmpvec
        }
      }
      curr_i = curr_i + 1
      next
    }

    # If the next value in the temporary vector is equal to the current value, then find the next value that is not equal to the current value
    start_i = curr_i
    while (tmpvec[curr_i] == curr_tmpvec) { curr_i = curr_i + 1 }

    # If the next value that is not equal to the current value is greater than the current value, then there is no peak near to the current value.
    # Move on to the next value that is greater than the current value.
    if (tmpvec[curr_i] > curr_tmpvec) { next }

    # If the next value that is not equal to the current value is less than the current value, then there is a peak near to the current value.
    # Count the peak and update the highest peak if necessary. Move on to the next value that is less than the current value.
    n_peaks = n_peaks + 1
    if (n_peaks == 1) {
      highest_peak_x = 0.5*(xvec[start_i - 1] + xvec[curr_i - 2])
      highest_peak_y = curr_tmpvec
    } else {
      if (curr_tmpvec > highest_peak_y) {
        highest_peak_x = 0.5*(xvec[start_i - 1] + xvec[curr_i - 2])
        highest_peak_y = curr_tmpvec
      }
    }
  }
}

# Return the output parameter values
return(list(n_peaks = n_peaks, highest_peak_x = highest_peak_x, highest_peak_y = highest_peak_y,
            n_troughs = n_troughs, lowest_trough_x = lowest_trough_x, lowest_trough_y = lowest_trough_y))
}


################################################################################################################################################
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

# Find the number of peaks in the first run of positive numbers in the "mu" curve while also finding the highest peak and its corresponding density
info_peaks_troughs = get_npeaks_ntroughs(reconstructed_model_fit$V2[index_mu_first_pos:index_mu_last_pos],
                                         reconstructed_model_fit$mu[index_mu_first_pos:index_mu_last_pos], -1.0)



#### ABOVE FULLY READ AND TESTED



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




a = 1

# Return the ....
return(a)
}
