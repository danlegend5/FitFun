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
tmpvec = rep_len(edge_value, N_p2)
tmpvec[2:(N_p2 - 1)] = yvec

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

# Return the number of peaks and troughs, and the coordinates of the highest peak and the lowest trough
return(list(n_peaks = n_peaks, highest_peak_x = highest_peak_x, highest_peak_y = highest_peak_y,
            n_troughs = n_troughs, lowest_trough_x = lowest_trough_x, lowest_trough_y = lowest_trough_y))
}


################################################################################################################################################
get_curve_properties = function(xvec, yvec, ind_hi, k_jam) {

# Description: For a fitted model component (e.g. for "sigma", "nu", or "tau") that has been reconstructed on a regular grid of density ranging from
#              zero to some positive value, this function computes approximate values for some useful properties of the curve in this range.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Set the default output parameter values
curve_properties = list(y_0 = NA, dydx_0 = NA, x_ymax = NA, y_max = NA, x_ymin = NA, y_min = NA, n_peaks = NA, n_troughs = NA, y_kjam = NA, dydx_kjam = NA)

# Compute some useful quantities
ngrid = length(xvec)
grid_density_step = xvec[2]

# Determine y at zero density
curve_properties$y_0 = yvec[1]

# Estimate the gradient of y (with respect to density) at zero density
curve_properties$dydx_0 = (yvec[2] - yvec[1])/grid_density_step

# If all of the values in the "mu" curve are non-positive, then return the computed properties of the y curve so far
if (is.na(ind_hi)) {
  return(curve_properties)
}

# If the first run of positive numbers in the "mu" curve ends before the end of the "mu" curve
if (ind_hi < ngrid) {

  # Estimate the gradient of y (with respect to density) at jam density
  ind_hi_p1 = ind_hi + 1
  curve_properties$dydx_kjam = (yvec[ind_hi_p1] - yvec[ind_hi])/grid_density_step

  # Estimate y at jam density
  curve_properties$y_kjam = yvec[ind_hi] + curve_properties$dydx_kjam*(k_jam - xvec[ind_hi])

  # Create a temporary vector containing the y values from zero to jam density
  tmp_yvec = yvec[1:ind_hi_p1]
  tmp_yvec[ind_hi_p1] = curve_properties$y_kjam

  # If all of the values in the y curve from zero to jam density are the same
  if (all(tmp_yvec == tmp_yvec[1])) {

    # Record that there are no peaks or troughs in the y curve from zero to jam density
    curve_properties$y_max = tmp_yvec[1]
    curve_properties$y_min = tmp_yvec[1]
    curve_properties$n_peaks = 0
    curve_properties$n_troughs = 0

  # If not all of the values in the y curve from zero to jam density are the same
  } else {

    # For the y curve from zero to jam density, find the number of peaks and troughs while also finding the coordinates of the highest peak and the
    # lowest trough
    tmp_xvec = xvec[1:ind_hi_p1]
    tmp_xvec[ind_hi_p1] = k_jam
    min_tmp_yvec = min(tmp_yvec)
    if (min_tmp_yvec >= 0.0) {
      edge_value = -1.0
    } else {
      edge_value = 1.1*min_tmp_yvec
    }
    info_peaks_troughs = get_npeaks_ntroughs(tmp_xvec, tmp_yvec, edge_value)
    curve_properties$x_ymax = info_peaks_troughs$highest_peak_x
    curve_properties$y_max = info_peaks_troughs$highest_peak_y
    curve_properties$n_peaks = info_peaks_troughs$n_peaks
    if (curve_properties$y_max <= 0.0) {
      edge_value = 1.0
    } else {
      edge_value = 1.1*curve_properties$y_max
    }
    info_peaks_troughs = get_npeaks_ntroughs(tmp_xvec, tmp_yvec, edge_value)
    curve_properties$x_ymin = info_peaks_troughs$lowest_trough_x
    curve_properties$y_min = info_peaks_troughs$lowest_trough_y
    curve_properties$n_troughs = info_peaks_troughs$n_troughs
  }

# If the first run of positive numbers in the "mu" curve does not end before the end of the "mu" curve
} else {

  # If all of the values in the y curve are the same
  if (all(yvec == yvec[1])) {

    # Record that there are no peaks or troughs in the y curve
    curve_properties$y_max = yvec[1]
    curve_properties$y_min = yvec[1]
    curve_properties$n_peaks = 0
    curve_properties$n_troughs = 0

  # If not all of the values in the y curve are the same
  } else {

    # Find the number of peaks and troughs in the y curve while also finding the coordinates of the highest peak and the lowest trough
    min_yvec = min(yvec)
    if (min_yvec >= 0.0) {
      edge_value = -1.0
    } else {
      edge_value = 1.1*min_yvec
    }
    info_peaks_troughs = get_npeaks_ntroughs(xvec, yvec, edge_value)
    curve_properties$x_ymax = info_peaks_troughs$highest_peak_x
    curve_properties$y_max = info_peaks_troughs$highest_peak_y
    curve_properties$n_peaks = info_peaks_troughs$n_peaks
    if (curve_properties$y_max <= 0.0) {
      edge_value = 1.0
    } else {
      edge_value = 1.1*curve_properties$y_max
    }
    info_peaks_troughs = get_npeaks_ntroughs(xvec, yvec, edge_value)
    curve_properties$x_ymin = info_peaks_troughs$lowest_trough_x
    curve_properties$y_min = info_peaks_troughs$lowest_trough_y
    curve_properties$n_troughs = info_peaks_troughs$n_troughs
  }
}

# Return the computed properties of the curve
return(curve_properties)
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
curve_properties_for_mu = list(q_0 = NA, v_ff = NA, dvdk_0 = NA, k_crit = NA, k_vmax = NA, q_cap = NA, v_max = NA, n_peaks = NA,
                               k_jam = NA, v_bw = NA, dvdk_kjam = NA, ind_first_pos = NA, ind_last_pos = NA)

# Compute some useful quantities
ngrid = nrow(reconstructed_model_fit)
grid_density_step = reconstructed_model_fit$V2[2]

# If the fitted model component for "mu" corresponds to flow
if (fd_type == 'Flow.Density') {

  # Determine the flow at zero density
  curve_properties_for_mu$q_0 = reconstructed_model_fit$mu[1]

  # Estimate the free-flow speed
  curve_properties_for_mu$v_ff = reconstructed_model_fit$mu[2]/grid_density_step

# If the fitted model component for "mu" corresponds to speed
} else if (fd_type == 'Speed.Density') {

  # Determine the free-flow speed
  curve_properties_for_mu$v_ff = reconstructed_model_fit$mu[1]

  # Estimate the gradient of the speed (with respect to density) at zero density
  curve_properties_for_mu$dvdk_0 = (reconstructed_model_fit$mu[2] - reconstructed_model_fit$mu[1])/grid_density_step
}

# Determine the grid index where the "mu" curve first becomes positive (if at all), and then determine the grid index at one step before the "mu"
# curve becomes non-positive again (if at all)
for (i in 1:ngrid) {
  if (reconstructed_model_fit$mu[i] <= 0.0) { next }
  curve_properties_for_mu$ind_first_pos = i
  curve_properties_for_mu$ind_last_pos = ngrid
  for (j in curve_properties_for_mu$ind_first_pos:ngrid) {
    if (reconstructed_model_fit$mu[j] > 0.0) { next }
    curve_properties_for_mu$ind_last_pos = j - 1
    break
  }
  break
}

# If all of the values in the "mu" curve are non-positive, then return the computed properties of the "mu" curve so far
if (is.na(curve_properties_for_mu$ind_first_pos)) {
  return(curve_properties_for_mu)
}

# Find the number of peaks in the first run of positive numbers in the "mu" curve while also finding the highest peak and its corresponding density
info_peaks_troughs = get_npeaks_ntroughs(reconstructed_model_fit$V2[curve_properties_for_mu$ind_first_pos:curve_properties_for_mu$ind_last_pos],
                                         reconstructed_model_fit$mu[curve_properties_for_mu$ind_first_pos:curve_properties_for_mu$ind_last_pos], -1.0)

# If the fitted model component for "mu" corresponds to flow
if (fd_type == 'Flow.Density') {

  # Record the critical density and the capacity
  curve_properties_for_mu$k_crit = info_peaks_troughs$highest_peak_x
  curve_properties_for_mu$q_cap = info_peaks_troughs$highest_peak_y

# If the fitted model component for "mu" corresponds to speed
} else if (fd_type == 'Speed.Density') {

  # Record the maximum speed and the corresponding density
  curve_properties_for_mu$k_vmax = info_peaks_troughs$highest_peak_x
  curve_properties_for_mu$v_max = info_peaks_troughs$highest_peak_y
}

# Record the number of peaks
curve_properties_for_mu$n_peaks = info_peaks_troughs$n_peaks

# If the first run of positive numbers in the "mu" curve ends before the end of the "mu" curve
if (curve_properties_for_mu$ind_last_pos < ngrid) {

  # Estimate the gradient (with respect to density) at the end of the first run of positive numbers in the "mu" curve
  grad = (reconstructed_model_fit$mu[curve_properties_for_mu$ind_last_pos + 1] - reconstructed_model_fit$mu[curve_properties_for_mu$ind_last_pos])/grid_density_step

  # Estimate the jam density
  curve_properties_for_mu$k_jam = reconstructed_model_fit$V2[curve_properties_for_mu$ind_last_pos] - (reconstructed_model_fit$mu[curve_properties_for_mu$ind_last_pos]/grad)

  # If the fitted model component for "mu" corresponds to flow
  if (fd_type == 'Flow.Density') {

    # Estimate the back-propagating wave speed at jam density
    curve_properties_for_mu$v_bw = -grad

    # Estimate the gradient of the speed (with respect to density) at jam density
    curve_properties_for_mu$dvdk_kjam = grad/curve_properties_for_mu$k_jam

  # If the fitted model component for "mu" corresponds to speed
  } else if (fd_type == 'Speed.Density') {

    # Estimate the back-propagating wave speed at jam density
    curve_properties_for_mu$v_bw = -grad*curve_properties_for_mu$k_jam

    # Estimate the gradient of the speed (with respect to density) at jam density
    curve_properties_for_mu$dvdk_kjam = grad
  }
}

# Return the computed properties of the "mu" curve
return(curve_properties_for_mu)
}


################################################################################################################################################
get_curve_properties_for_sigma = function(reconstructed_model_fit, curve_properties_for_mu) {

# Description: For a fitted model component for "sigma" that has been reconstructed on a regular grid of density ranging from zero to some positive
#              value, this function computes approximate values for some useful properties of the curve in this range.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Compute approximate values for some useful properties of the "sigma" curve
curve_properties = get_curve_properties(reconstructed_model_fit$V2, reconstructed_model_fit$sigma, curve_properties_for_mu$ind_last_pos,
                                        curve_properties_for_mu$k_jam)

# Return the computed properties of the "sigma" curve
return(list(sigma_0 = curve_properties$y_0, dsigmadk_0 = curve_properties$dydx_0, k_sigmamax = curve_properties$x_ymax,
            sigma_max = curve_properties$y_max, k_sigmamin = curve_properties$x_ymin, sigma_min = curve_properties$y_min,
            n_peaks = curve_properties$n_peaks, n_troughs = curve_properties$n_troughs, sigma_kjam = curve_properties$y_kjam,
            dsigmadk_kjam = curve_properties$dydx_kjam))
}


################################################################################################################################################
get_curve_properties_for_nu = function(reconstructed_model_fit, curve_properties_for_mu) {

# Description: For a fitted model component for "nu" that has been reconstructed on a regular grid of density ranging from zero to some positive
#              value, this function computes approximate values for some useful properties of the curve in this range.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Compute approximate values for some useful properties of the "nu" curve
curve_properties = get_curve_properties(reconstructed_model_fit$V2, reconstructed_model_fit$nu, curve_properties_for_mu$ind_last_pos,
                                        curve_properties_for_mu$k_jam)

# Return the computed properties of the "nu" curve
return(list(nu_0 = curve_properties$y_0, dnudk_0 = curve_properties$dydx_0, k_numax = curve_properties$x_ymax,
            nu_max = curve_properties$y_max, k_numin = curve_properties$x_ymin, nu_min = curve_properties$y_min,
            n_peaks = curve_properties$n_peaks, n_troughs = curve_properties$n_troughs, nu_kjam = curve_properties$y_kjam,
            dnudk_kjam = curve_properties$dydx_kjam))
}


################################################################################################################################################
get_curve_properties_for_tau = function(reconstructed_model_fit, curve_properties_for_mu) {

# Description: For a fitted model component for "tau" that has been reconstructed on a regular grid of density ranging from zero to some positive
#              value, this function computes approximate values for some useful properties of the curve in this range.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Compute approximate values for some useful properties of the "tau" curve
curve_properties = get_curve_properties(reconstructed_model_fit$V2, reconstructed_model_fit$tau, curve_properties_for_mu$ind_last_pos,
                                        curve_properties_for_mu$k_jam)

# Return the computed properties of the "tau" curve
return(list(tau_0 = curve_properties$y_0, dtaudk_0 = curve_properties$dydx_0, k_taumax = curve_properties$x_ymax,
            tau_max = curve_properties$y_max, k_taumin = curve_properties$x_ymin, tau_min = curve_properties$y_min,
            n_peaks = curve_properties$n_peaks, n_troughs = curve_properties$n_troughs, tau_kjam = curve_properties$y_kjam,
            dtaudk_kjam = curve_properties$dydx_kjam))
}


################################################################################################################################################
write_fit_summary = function(output_file, fd_type, ndata, data_min_density, data_max_density, data_min_y, data_max_y,
                             npar_mu, npar_sigma, npar_nu, npar_tau, npar_all, gdev, aic, bic,
                             q_0, v_ff, dvdk_0, k_crit, k_vmax, q_cap, v_max, k_jam, v_bw, dvdk_kjam,
                             curve_properties_for_mu_over_data_range, curve_properties_for_sigma_over_data_range,
                             curve_properties_for_nu_over_data_range, curve_properties_for_tau_over_data_range,
                             curve_properties_for_mu_over_full_range, curve_properties_for_sigma_over_full_range,
                             curve_properties_for_nu_over_full_range, curve_properties_for_tau_over_full_range) {

# Description: This function writes out a fit summary to the file "output_file".
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Write out the fit summary file "output_file"
cat('######################################################################################################################\n', file = output_file)
cat('# DATA PROPERTIES\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
if (fd_type == 'Flow.Density') {
  cat(ndata, '          # No. of flow-density measurement pairs (Ndat)\n', file = output_file, append = TRUE)
} else if (fd_type == 'Speed.Density') {
  cat(ndata, '          # No. of speed-density measurement pairs (Ndat)\n', file = output_file, append = TRUE)
}
cat(data_min_density, '          # Minimum density in the data\n', file = output_file, append = TRUE)
cat(data_max_density, '          # Maximum density in the data\n', file = output_file, append = TRUE)
if (fd_type == 'Flow.Density') {
  cat(data_min_y, '          # Minimum flow in the data\n', file = output_file, append = TRUE)
  cat(data_max_y, '          # Maximum flow in the data\n', file = output_file, append = TRUE)
} else if (fd_type == 'Speed.Density') {
  cat(data_min_y, '          # Minimum speed in the data\n', file = output_file, append = TRUE)
  cat(data_max_y, '          # Maximum speed in the data\n', file = output_file, append = TRUE)
}
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# MODEL PARAMETER COUNTS\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(npar_mu, '          # No. of free parameters for "mu"\n', file = output_file, append = TRUE)
cat(npar_sigma, '          # No. of free parameters for "sigma"\n', file = output_file, append = TRUE)
cat(npar_nu, '          # No. of free parameters for "nu"\n', file = output_file, append = TRUE)
cat(npar_tau, '          # No. of free parameters for "tau"\n', file = output_file, append = TRUE)
cat(npar_all, '          # Total no. of free parameters (Npar)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# FIT QUALITY\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(gdev, '          # Global deviance (-2 ln L)\n', file = output_file, append = TRUE)
cat(aic, '          # AIC (-2 ln L + 2 Npar)\n', file = output_file, append = TRUE)
cat(bic, '          # BIC (-2 ln L + Npar ln Ndat)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# FITTED PHYSICAL PARAMETERS\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(q_0, '          # Flow at zero density\n', file = output_file, append = TRUE)
cat(v_ff, '          # Free-flow speed (i.e. the speed as density tends to zero)\n', file = output_file, append = TRUE)
cat(dvdk_0, '          # Gradient of the speed (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(k_crit, '          # Critical density (i.e. the density at maximum flow)\n', file = output_file, append = TRUE)
cat(k_vmax, '          # Density at maximum speed\n', file = output_file, append = TRUE)
cat(q_cap, '          # Capacity (i.e. the maximum flow)\n', file = output_file, append = TRUE)
cat(v_max, '          # Maximum speed\n', file = output_file, append = TRUE)
cat(k_jam, '          # Jam density (i.e. the density at zero speed)\n', file = output_file, append = TRUE)
cat(v_bw, '          # Back-propagating wave speed at jam density\n', file = output_file, append = TRUE)
cat(dvdk_kjam, '          # Gradient of the speed (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "MU" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE FIRST RUN OF POSITIVE NUMBERS IN THE "MU" CURVE\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$q_0, '          # Flow at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$v_ff, '          # Free-flow speed (i.e. the speed as density tends to zero)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$dvdk_0, '          # Gradient of the speed (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$k_crit, '          # *Critical density (i.e. the density at maximum flow)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$k_vmax, '          # *Density at maximum speed\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$q_cap, '          # *Capacity (i.e. the maximum flow)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$v_max, '          # *Maximum speed\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$k_jam, '          # *Jam density (i.e. the density at zero speed)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$v_bw, '          # *Back-propagating wave speed at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_data_range$dvdk_kjam, '          # *Gradient of the speed (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "SIGMA" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$sigma_0, '          # Sigma at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$dsigmadk_0, '          # Gradient of the sigma (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$k_sigmamax, '          # *Density at maximum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$sigma_max, '          # *Maximum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$k_sigmamin, '          # *Density at minimum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$sigma_min, '          # *Minimum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$n_troughs, '          # *No. of troughs\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$sigma_kjam, '          # Sigma at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_data_range$dsigmadk_kjam, '          # Gradient of the sigma (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "NU" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$nu_0, '          # Nu at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$dnudk_0, '          # Gradient of the nu (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$k_numax, '          # *Density at maximum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$nu_max, '          # *Maximum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$k_numin, '          # *Density at minimum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$nu_min, '          # *Minimum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$n_troughs, '          # *No. of troughs\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$nu_kjam, '          # Nu at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_data_range$dnudk_kjam, '          # Gradient of the nu (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "TAU" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$tau_0, '          # Tau at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$dtaudk_0, '          # Gradient of the tau (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$k_taumax, '          # *Density at maximum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$tau_max, '          # *Maximum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$k_taumin, '          # *Density at minimum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$tau_min, '          # *Minimum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$n_troughs, '          # *No. of troughs\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$tau_kjam, '          # Tau at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_data_range$dtaudk_kjam, '          # Gradient of the tau (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "MU" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE FIRST RUN OF POSITIVE NUMBERS IN THE "MU" CURVE\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$q_0, '          # Flow at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$v_ff, '          # Free-flow speed (i.e. the speed as density tends to zero)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$dvdk_0, '          # Gradient of the speed (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$k_crit, '          # *Critical density (i.e. the density at maximum flow)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$k_vmax, '          # *Density at maximum speed\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$q_cap, '          # *Capacity (i.e. the maximum flow)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$v_max, '          # *Maximum speed\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$k_jam, '          # *Jam density (i.e. the density at zero speed)\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$v_bw, '          # *Back-propagating wave speed at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_mu_over_full_range$dvdk_kjam, '          # *Gradient of the speed (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "SIGMA" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$sigma_0, '          # Sigma at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$dsigmadk_0, '          # Gradient of the sigma (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$k_sigmamax, '          # *Density at maximum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$sigma_max, '          # *Maximum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$k_sigmamin, '          # *Density at minimum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$sigma_min, '          # *Minimum sigma\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$n_troughs, '          # *No. of troughs\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$sigma_kjam, '          # Sigma at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_sigma_over_full_range$dsigmadk_kjam, '          # Gradient of the sigma (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "NU" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$nu_0, '          # Nu at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$dnudk_0, '          # Gradient of the nu (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$k_numax, '          # *Density at maximum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$nu_max, '          # *Maximum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$k_numin, '          # *Density at minimum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$nu_min, '          # *Minimum nu\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$n_troughs, '          # *No. of troughs\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$nu_kjam, '          # Nu at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_nu_over_full_range$dnudk_kjam, '          # Gradient of the nu (with respect to density) at jam density\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat('# ESTIMATED PROPERTIES OF THE FITTED "TAU" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n', file = output_file, append = TRUE)
cat('# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n', file = output_file, append = TRUE)
cat('######################################################################################################################\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$tau_0, '          # Tau at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$dtaudk_0, '          # Gradient of the tau (with respect to density) at zero density\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$k_taumax, '          # *Density at maximum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$tau_max, '          # *Maximum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$k_taumin, '          # *Density at minimum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$tau_min, '          # *Minimum tau\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$n_peaks, '          # *No. of peaks\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$n_troughs, '          # *No. of troughs\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$tau_kjam, '          # Tau at jam density\n', file = output_file, append = TRUE)
cat(curve_properties_for_tau_over_full_range$dtaudk_kjam, '          # Gradient of the tau (with respect to density) at jam density\n', file = output_file, append = TRUE)
}
