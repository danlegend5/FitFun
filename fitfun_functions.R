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
  # previous value in the temporary vector is equal to the current value).
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

# Description: For a fitted model component (e.g. for "sigma", "nu", or "tau") that has been reconstructed on a regular grid of density ranging
#              from zero to some positive value, this function computes approximate values for some useful properties of the curve in this range.
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

# Determine the value of the curve at zero density
curve_properties$y_0 = yvec[1]

# Estimate the gradient of the curve (with respect to density) at zero density
curve_properties$dydx_0 = (yvec[2] - yvec[1])/grid_density_step

# If all of the values in the curve are non-positive, then return the computed properties of the curve so far
if (is.na(ind_hi)) {
  return(curve_properties)
}

# If the first run of positive numbers in the curve ends before the end of the curve
if (ind_hi < ngrid) {

  # Estimate the gradient of the curve (with respect to density) at jam density
  ind_hi_p1 = ind_hi + 1
  curve_properties$dydx_kjam = (yvec[ind_hi_p1] - yvec[ind_hi])/grid_density_step

  # Estimate the value of the curve at jam density
  curve_properties$y_kjam = yvec[ind_hi] + curve_properties$dydx_kjam*(k_jam - xvec[ind_hi])

  # Create a temporary vector containing the values of the curve from zero to jam density
  tmp_yvec = yvec[1:ind_hi_p1]
  tmp_yvec[ind_hi_p1] = curve_properties$y_kjam

  # If all of the values of the curve from zero to jam density are the same
  if (all(tmp_yvec == tmp_yvec[1])) {

    # Record that there are no peaks or troughs in the curve from zero to jam density
    curve_properties$y_max = tmp_yvec[1]
    curve_properties$y_min = tmp_yvec[1]
    curve_properties$n_peaks = 0
    curve_properties$n_troughs = 0

  # If not all of the values of the curve from zero to jam density are the same
  } else {

    # For the curve from zero to jam density, find the number of peaks and troughs while also finding the coordinates of the highest peak and the
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

# If the first run of positive numbers in the curve does not end before the end of the curve
} else {

  # If all of the values of the curve are the same
  if (all(yvec == yvec[1])) {

    # Record that there are no peaks or troughs in the curve
    curve_properties$y_max = yvec[1]
    curve_properties$y_min = yvec[1]
    curve_properties$n_peaks = 0
    curve_properties$n_troughs = 0

  # If not all of the values of the curve are the same
  } else {

    # Find the number of peaks and troughs in the curve while also finding the coordinates of the highest peak and the lowest trough
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

  # Estimate the gradient of the speed (with respect to density) at zero density
  if (ngrid > 2) {
    curve_properties_for_mu$dvdk_0 = ((reconstructed_model_fit$mu[3]/reconstructed_model_fit$V2[3]) - curve_properties_for_mu$v_ff)/grid_density_step
  }

# If the fitted model component for "mu" corresponds to speed
} else if (fd_type == 'Speed.Density') {

  # Determine the flow at zero density
  curve_properties_for_mu$q_0 = 0.0

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

# Record the number of peaks
curve_properties_for_mu$n_peaks = info_peaks_troughs$n_peaks

# If the fitted model component for "mu" corresponds to flow
if (fd_type == 'Flow.Density') {

  # Record the critical density and the capacity
  curve_properties_for_mu$k_crit = info_peaks_troughs$highest_peak_x
  curve_properties_for_mu$q_cap = info_peaks_troughs$highest_peak_y

  # Construct a speed curve, and estimate the maximum speed and the corresponding density
  if (curve_properties_for_mu$ind_last_pos > 1) {
    ind_lo = max(2, curve_properties_for_mu$ind_first_pos)
    tmp_density_vec = reconstructed_model_fit$V2[ind_lo:curve_properties_for_mu$ind_last_pos]
    speed_curve = reconstructed_model_fit$mu[ind_lo:curve_properties_for_mu$ind_last_pos]/tmp_density_vec
    info_peaks_troughs = get_npeaks_ntroughs(tmp_density_vec, speed_curve, -1.0)
    curve_properties_for_mu$k_vmax = info_peaks_troughs$highest_peak_x
    curve_properties_for_mu$v_max = info_peaks_troughs$highest_peak_y
  }

# If the fitted model component for "mu" corresponds to speed
} else if (fd_type == 'Speed.Density') {

  # Record the maximum speed and the corresponding density
  curve_properties_for_mu$k_vmax = info_peaks_troughs$highest_peak_x
  curve_properties_for_mu$v_max = info_peaks_troughs$highest_peak_y

  # Construct a flow curve, and estimate the critical density and the capacity
  tmp_density_vec = reconstructed_model_fit$V2[curve_properties_for_mu$ind_first_pos:curve_properties_for_mu$ind_last_pos]
  flow_curve = reconstructed_model_fit$mu[curve_properties_for_mu$ind_first_pos:curve_properties_for_mu$ind_last_pos]*tmp_density_vec
  info_peaks_troughs = get_npeaks_ntroughs(tmp_density_vec, flow_curve, -1.0)
  curve_properties_for_mu$k_crit = info_peaks_troughs$highest_peak_x
  curve_properties_for_mu$q_cap = info_peaks_troughs$highest_peak_y
}

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
calculate_normalised_quantile_residuals = function(cprobs_lower, cprobs_upper) {

# Description: This function converts lower cumulative probabilities "cprobs_lower" and upper cumulative probabilities "cprobs_upper" to standard
#              Normal quantiles using the inverse cumulative distribution function (CDF) of the standard Normal distribution. If these cumulative
#              probabilities correspond to the CDF values of a fitted GAMLSS model at the values of a set of response observations, then the
#              output standard Normal quantiles are referred to as normalised quantile residuals. This function has been implemented to provide
#              normalised quantile residuals with less occurrences of "-Inf" and "Inf" values than those provided by the "gamlss()" function.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)


# Calculate the normalised quantile residuals
ndata = length(cprobs_lower)
selection = cprobs_lower < 0.5
nselection = sum(selection)
if (nselection == ndata) {
  nqr = qNO(cprobs_lower)
} else {
  if (nselection == 0) {
    nqr = qNO(cprobs_upper, lower.tail = FALSE)
  } else {
    nqr = double(length = ndata)
    nqr[selection] = qNO(cprobs_lower[selection])
    selection = !selection
    nqr[selection] = qNO(cprobs_upper[selection], lower.tail = FALSE)
  }
}

# Return the normalised quantile residuals
return(nqr)
}


################################################################################################################################################
calculate_distributional_measures_for_SN2 = function(mu, sigma, nu) {

# Description: Given the location, scale, and skewness parameters "mu", "sigma", and "nu", respectively, of a Skew Normal Type II distribution,
#              this function calculates the corresponding mean, median, mode, standard deviation, moment skewness, and moment excess kurtosis of
#              the distribution. These distributional measures are calculated using the formulae given in Rigby et al. (2019, CRC Press LLC, New
#              York).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)


# Calculate the mean of the distribution
fac = sqrt(2.0/pi)
inv_nu = 1.0/nu
ez = fac*(nu - inv_nu)
mean = mu + (sigma*ez)

# Calculate the median of the distribution
nu2 = nu*nu
median = (1.0 + nu2)/4.0
selection = nu > 1.0
nselection = sum(selection)
if (nselection > 0) {
  tmpvec = nu2[selection]
  median[selection] = ((3.0*tmpvec) - 1.0)/(4.0*tmpvec)
}
median = (sigma/nu)*qNO(median)
if (nselection > 0) {
  median[selection] = median[selection]*tmpvec
}
median = median + mu

# Calculate the mode of the distribution
mode = mu

# Calculate the standard deviation of the distribution
inv_nu2 = inv_nu*inv_nu
ez2 = ez*ez
varz = nu2 + inv_nu2 - 1.0 - ez2
stddev = sigma*sqrt(varz)

# Calculate the moment skewness of the distribution
nu_p_inv_nu = nu + inv_nu
nu4 = nu2*nu2
inv_nu4 = inv_nu2*inv_nu2
tmpvec1 = (2.0*fac)*((nu4 - inv_nu4)/nu_p_inv_nu)
skewness = (tmpvec1 - (3.0*varz*ez) - (ez2*ez))/(varz^1.5)

# Calculate the moment excess kurtosis of the distribution
tmpvec2 = 3.0*((nu4*nu) + (inv_nu4*inv_nu))/nu_p_inv_nu
excess_kurtosis = ((tmpvec2 - (4.0*tmpvec1*ez) + (6.0*varz*ez2) + (3.0*ez2*ez2))/(varz*varz)) - 3.0

# Return the distributional measures
return(data.table(mean = mean, median = median, mode = mode, standard_deviation = stddev, moment_skewness = skewness, moment_excess_kurtosis = excess_kurtosis))
}


################################################################################################################################################
calculate_distributional_measures_for_SEP3 = function(mu, sigma, nu, tau) {

# Description: Given the location, scale, skewness, and kurtosis parameters "mu", "sigma", "nu", and "tau", respectively, of a Skew Exponential
#              Power Type III distribution, this function calculates the corresponding mean, mode, standard deviation, moment skewness, and moment
#              excess kurtosis of the distribution. These distributional measures are calculated using the formulae given in Rigby et al. (2019, CRC
#              Press LLC, New York). This reference does not supply an analytical formula for the median which could be because it does not exist.
#              Consequently, this function does not calculate the median.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)


# Calculate the mean of the distribution
inv_nu = 1.0/nu
inv_tau = 1.0/tau
two_power_inv_tau = 2.0^inv_tau
gamma_inv_tau = gamma(inv_tau)
ez = two_power_inv_tau*gamma(2.0*inv_tau)*(nu - inv_nu)/gamma_inv_tau
mean = mu + (sigma*ez)

# Calculate the mode of the distribution
mode = mu

# Calculate the standard deviation of the distribution
two_power_inv_tau_squared = two_power_inv_tau*two_power_inv_tau
nu2 = nu*nu
inv_nu2 = inv_nu*inv_nu
ez2 = ez*ez
varz = (two_power_inv_tau_squared*gamma(3.0*inv_tau)*(nu2 + inv_nu2 - 1.0)/gamma_inv_tau) - ez2
stddev = sigma*sqrt(varz)

# Calculate the moment skewness of the distribution
nu4 = nu2*nu2
inv_nu4 = inv_nu2*inv_nu2
denom = gamma_inv_tau*(nu + inv_nu)
tmpvec1 = two_power_inv_tau_squared*two_power_inv_tau*gamma(4.0*inv_tau)*(nu4 - inv_nu4)/denom
skewness = (tmpvec1 - (3.0*varz*ez) - (ez2*ez))/(varz^1.5)

# Calculate the moment excess kurtosis of the distribution
tmpvec2 = two_power_inv_tau_squared*two_power_inv_tau_squared*gamma(5.0*inv_tau)*((nu4*nu) + (inv_nu4*inv_nu))/denom
excess_kurtosis = ((tmpvec2 - (4.0*tmpvec1*ez) + (6.0*varz*ez2) + (3.0*ez2*ez2))/(varz*varz)) - 3.0

# Return the distributional measures
return(data.table(mean = mean, mode = mode, standard_deviation = stddev, moment_skewness = skewness, moment_excess_kurtosis = excess_kurtosis))
}


################################################################################################################################################
fix_out_of_data_curves = function(xvec, yvec, min_x, max_x) {

# Description: For a curve in the xy-plane that is represented by a discrete set of N points with (x,y) pairs "xvec" and "yvec", this function does
#              the following. Firstly, it determines the y-value "y_lower" corresponding to the minimum x-value that is greater than or equal to
#              "min_x", and then, for all x-values less than "min_x", it sets the y-values in the curve to "y_lower". Secondly, it determines the
#              y-value "y_upper" corresponding to the maximum x-value that is less than or equal to "max_x", and then, for all x-values greater than
#              "max_x", it sets the y-values in the curve to "y_upper".
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Determine the set of points with x-values in the range "min_x" to "max_x"
selection = (xvec >= min_x) & (xvec <= max_x)

# If there are no points with x-values in the range "min_x" to "max_x", then do nothing and return
if (sum(selection) == 0) { return(yvec) }

# Extract the y-values at the end points of the selection
tmpvec = yvec[selection]
y_lower = tmpvec[1]
y_upper = tmpvec[length(tmpvec)]

# For all x-values less than "min_x", set the y-values in the curve to "y_lower"
selection = xvec < min_x
if (sum(selection) > 0) { yvec[selection] = y_lower }

# For all x-values greater than "max_x", set the y-values in the curve to "y_upper"
selection = xvec > max_x
if (sum(selection) > 0) { yvec[selection] = y_upper }

# Return the modified y-values
return(yvec)
}


################################################################################################################################################
compute_slotted_acf = function(tvec, yvec, tlag_bin_size, tlag_nbins) {

# Description: This function computes the slotted auto-correlation function (slotted ACF) for a set of time series data "yvec" observed at irregular
#              times "tvec". The time series data must be sorted such that the values in "tvec" are in ascending order. The value of the slotted ACF
#              for a specific time lag bin is computed by averaging the cross products of sample pairs whose time differences fall in the given bin
#              (Edelson & Krolik, 1988, Astrophysical Journal, 333, 646). As part of the computation, this function employs the algorithm for
#              calculating the running mean and variance invented by Welford (1962, Technometrics, 4, 419). The time lag bin size "tlag_bin_size" in
#              combination with the number of time lag bins "tlag_nbins" specifies the maximum time difference to consider in the computation of the
#              slotted ACF.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Determine the time lag bin limits
tlag_upper_bin_limits = tlag_bin_size*seq(from = 0.5, to = tlag_nbins - 0.5, length.out = tlag_nbins)

# Set up some vectors and arrays
ndata = length(tvec)
npairs_per_bin = double(length = tlag_nbins)
obsflag_per_bin = array(0, dim = c(tlag_nbins, ndata))
acf = double(length = tlag_nbins)
acf_err = double(length = tlag_nbins)

# For each time series observation
max_tlag = tlag_upper_bin_limits[tlag_nbins]
for (curr_i in 1:ndata) {

  # Extract the current data values
  t_curr_i = tvec[curr_i]
  y_curr_i = yvec[curr_i]

  # Determine the cutoff time
  t_cutoff = t_curr_i + max_tlag

  # Step forwards through the time series observations from the current observation
  curr_j = curr_i
  curr_bin = 1
  while (curr_j <= ndata) {

    # If the time of the future observation is beyond the cutoff time, then break out of the loop
    t_curr_j = tvec[curr_j]
    if (t_curr_j >= t_cutoff) { break }

    # Calculate the time lag
    tlag = t_curr_j - t_curr_i

    # Find the corresponding time lag bin
    while (curr_bin <= tlag_nbins) {
      if (tlag < tlag_upper_bin_limits[curr_bin]) { break }
      curr_bin = curr_bin + 1
    }

    # Compute the product of the observations
    y_cross = y_curr_i*yvec[curr_j]

    # Update the relevant vectors and arrays
    new_val_npairs = npairs_per_bin[curr_bin] + 1.0
    npairs_per_bin[curr_bin] = new_val_npairs
    obsflag_per_bin[curr_bin, curr_i] = 1
    if (new_val_npairs == 1.0) {
      acf[curr_bin] = y_cross
    } else {
      prev_val_acf = acf[curr_bin]
      tmp_resid = y_cross - prev_val_acf
      new_val_acf = prev_val_acf + (tmp_resid/new_val_npairs)
      acf[curr_bin] = new_val_acf
      acf_err[curr_bin] = acf_err[curr_bin] + (tmp_resid*(y_cross - new_val_acf))
    }

    # Move on to the next future observation
    curr_j = curr_j + 1
  }
}

# Compute the uncertainties on the slotted ACF values
obsflag_per_bin = apply(obsflag_per_bin, 1, sum)
selection_good = obsflag_per_bin > 1
ngood = sum(selection_good)
if (ngood == tlag_nbins) {
  acf_err = sqrt(acf_err/((npairs_per_bin - 1.0)*(obsflag_per_bin - 1.0)))
} else if (ngood == 0) {
  acf[1:tlag_nbins] = 0.0
  acf_err[1:tlag_nbins] = -1.0
} else {
  acf_err[selection_good] = sqrt(acf_err[selection_good]/((npairs_per_bin[selection_good] - 1.0)*(obsflag_per_bin[selection_good] - 1.0)))
  selection_bad = obsflag_per_bin <= 1
  acf[selection_bad] = 0.0
  acf_err[selection_bad] = -1.0
}

# Return the slotted ACF
return(data.table(tlag_bin_mid = tlag_upper_bin_limits - 0.5*tlag_bin_size, tlag_bin_lo = tlag_upper_bin_limits - tlag_bin_size,
                  tlag_bin_hi = tlag_upper_bin_limits, acf = acf, acf_err = acf_err))
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
cat('######################################################################################################################\n',
    '# DATA PROPERTIES\n',
    '######################################################################################################################\n',
    file = output_file, sep = '')
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
cat('######################################################################################################################\n',
    '# MODEL PARAMETER COUNTS\n',
    '######################################################################################################################\n',
    npar_mu, '           # No. of free parameters for "mu"\n',
    npar_sigma, '           # No. of free parameters for "sigma"\n',
    npar_nu, '           # No. of free parameters for "nu"\n',
    npar_tau, '           # No. of free parameters for "tau"\n',
    npar_all, '           # Total no. of free parameters (Npar)\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# FIT QUALITY\n',
    '######################################################################################################################\n',
    format(gdev, nsmall = 4), '           # Global deviance [ -2*ln(Lmax) ]\n',
    format(aic, nsmall = 4), '           # AIC [ -2*ln(Lmax) + 2*Npar ]\n',
    format(bic, nsmall = 4), '           # BIC [ -2*ln(Lmax) + Npar*ln(Ndat) ]\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# FITTED PHYSICAL PARAMETERS\n',
    '######################################################################################################################\n',
    q_0, '           # Flow at zero density\n',
    v_ff, '           # Free-flow speed (i.e. the speed as density tends to zero)\n',
    dvdk_0, '           # Gradient of the speed (with respect to density) at zero density\n',
    k_crit, '           # Critical density (i.e. the density at maximum flow)\n',
    k_vmax, '           # Density at maximum speed\n',
    q_cap, '           # Capacity (i.e. the maximum flow)\n',
    v_max, '           # Maximum speed\n',
    k_jam, '           # Jam density (i.e. the density at zero speed)\n',
    v_bw, '           # Back-propagating wave speed at jam density\n',
    dvdk_kjam, '           # Gradient of the speed (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "MU" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n',
    '# * = ESTIMATED FOR THE FIRST RUN OF POSITIVE NUMBERS IN THE "MU" CURVE\n',
    '######################################################################################################################\n',
    curve_properties_for_mu_over_data_range$q_0, '           # Flow at zero density\n',
    curve_properties_for_mu_over_data_range$v_ff, '           # Free-flow speed (i.e. the speed as density tends to zero)\n',
    curve_properties_for_mu_over_data_range$dvdk_0, '           # Gradient of the speed (with respect to density) at zero density\n',
    curve_properties_for_mu_over_data_range$k_crit, '           # *Critical density (i.e. the density at maximum flow)\n',
    curve_properties_for_mu_over_data_range$k_vmax, '           # *Density at maximum speed\n',
    curve_properties_for_mu_over_data_range$q_cap, '           # *Capacity (i.e. the maximum flow)\n',
    curve_properties_for_mu_over_data_range$v_max, '           # *Maximum speed\n',
    curve_properties_for_mu_over_data_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_mu_over_data_range$k_jam, '           # *Jam density (i.e. the density at zero speed)\n',
    curve_properties_for_mu_over_data_range$v_bw, '           # *Back-propagating wave speed at jam density\n',
    curve_properties_for_mu_over_data_range$dvdk_kjam, '           # *Gradient of the speed (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "SIGMA" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n',
    '# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n',
    '######################################################################################################################\n',
    curve_properties_for_sigma_over_data_range$sigma_0, '           # Sigma at zero density\n',
    curve_properties_for_sigma_over_data_range$dsigmadk_0, '           # Gradient of the sigma (with respect to density) at zero density\n',
    curve_properties_for_sigma_over_data_range$k_sigmamax, '           # *Density at maximum sigma\n',
    curve_properties_for_sigma_over_data_range$sigma_max, '           # *Maximum sigma\n',
    curve_properties_for_sigma_over_data_range$k_sigmamin, '           # *Density at minimum sigma\n',
    curve_properties_for_sigma_over_data_range$sigma_min, '           # *Minimum sigma\n',
    curve_properties_for_sigma_over_data_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_sigma_over_data_range$n_troughs, '           # *No. of troughs\n',
    curve_properties_for_sigma_over_data_range$sigma_kjam, '           # Sigma at jam density\n',
    curve_properties_for_sigma_over_data_range$dsigmadk_kjam, '           # Gradient of the sigma (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "NU" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n',
    '# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n',
    '######################################################################################################################\n',
    curve_properties_for_nu_over_data_range$nu_0, '           # Nu at zero density\n',
    curve_properties_for_nu_over_data_range$dnudk_0, '           # Gradient of the nu (with respect to density) at zero density\n',
    curve_properties_for_nu_over_data_range$k_numax, '           # *Density at maximum nu\n',
    curve_properties_for_nu_over_data_range$nu_max, '           # *Maximum nu\n',
    curve_properties_for_nu_over_data_range$k_numin, '           # *Density at minimum nu\n',
    curve_properties_for_nu_over_data_range$nu_min, '           # *Minimum nu\n',
    curve_properties_for_nu_over_data_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_nu_over_data_range$n_troughs, '           # *No. of troughs\n',
    curve_properties_for_nu_over_data_range$nu_kjam, '           # Nu at jam density\n',
    curve_properties_for_nu_over_data_range$dnudk_kjam, '           # Gradient of the nu (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "TAU" MODEL OVER THE DENSITY RANGE FROM ZERO TO THE MAXIMUM DENSITY IN THE DATA\n',
    '# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n',
    '######################################################################################################################\n',
    curve_properties_for_tau_over_data_range$tau_0, '           # Tau at zero density\n',
    curve_properties_for_tau_over_data_range$dtaudk_0, '           # Gradient of the tau (with respect to density) at zero density\n',
    curve_properties_for_tau_over_data_range$k_taumax, '           # *Density at maximum tau\n',
    curve_properties_for_tau_over_data_range$tau_max, '           # *Maximum tau\n',
    curve_properties_for_tau_over_data_range$k_taumin, '           # *Density at minimum tau\n',
    curve_properties_for_tau_over_data_range$tau_min, '           # *Minimum tau\n',
    curve_properties_for_tau_over_data_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_tau_over_data_range$n_troughs, '           # *No. of troughs\n',
    curve_properties_for_tau_over_data_range$tau_kjam, '           # Tau at jam density\n',
    curve_properties_for_tau_over_data_range$dtaudk_kjam, '           # Gradient of the tau (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "MU" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n',
    '# * = ESTIMATED FOR THE FIRST RUN OF POSITIVE NUMBERS IN THE "MU" CURVE\n',
    '######################################################################################################################\n',
    curve_properties_for_mu_over_full_range$q_0, '           # Flow at zero density\n',
    curve_properties_for_mu_over_full_range$v_ff, '           # Free-flow speed (i.e. the speed as density tends to zero)\n',
    curve_properties_for_mu_over_full_range$dvdk_0, '           # Gradient of the speed (with respect to density) at zero density\n',
    curve_properties_for_mu_over_full_range$k_crit, '           # *Critical density (i.e. the density at maximum flow)\n',
    curve_properties_for_mu_over_full_range$k_vmax, '           # *Density at maximum speed\n',
    curve_properties_for_mu_over_full_range$q_cap, '           # *Capacity (i.e. the maximum flow)\n',
    curve_properties_for_mu_over_full_range$v_max, '           # *Maximum speed\n',
    curve_properties_for_mu_over_full_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_mu_over_full_range$k_jam, '           # *Jam density (i.e. the density at zero speed)\n',
    curve_properties_for_mu_over_full_range$v_bw, '           # *Back-propagating wave speed at jam density\n',
    curve_properties_for_mu_over_full_range$dvdk_kjam, '           # *Gradient of the speed (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "SIGMA" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n',
    '# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n',
    '######################################################################################################################\n',
    curve_properties_for_sigma_over_full_range$sigma_0, '           # Sigma at zero density\n',
    curve_properties_for_sigma_over_full_range$dsigmadk_0, '           # Gradient of the sigma (with respect to density) at zero density\n',
    curve_properties_for_sigma_over_full_range$k_sigmamax, '           # *Density at maximum sigma\n',
    curve_properties_for_sigma_over_full_range$sigma_max, '           # *Maximum sigma\n',
    curve_properties_for_sigma_over_full_range$k_sigmamin, '           # *Density at minimum sigma\n',
    curve_properties_for_sigma_over_full_range$sigma_min, '           # *Minimum sigma\n',
    curve_properties_for_sigma_over_full_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_sigma_over_full_range$n_troughs, '           # *No. of troughs\n',
    curve_properties_for_sigma_over_full_range$sigma_kjam, '           # Sigma at jam density\n',
    curve_properties_for_sigma_over_full_range$dsigmadk_kjam, '           # Gradient of the sigma (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "NU" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n',
    '# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n',
    '######################################################################################################################\n',
    curve_properties_for_nu_over_full_range$nu_0, '           # Nu at zero density\n',
    curve_properties_for_nu_over_full_range$dnudk_0, '           # Gradient of the nu (with respect to density) at zero density\n',
    curve_properties_for_nu_over_full_range$k_numax, '           # *Density at maximum nu\n',
    curve_properties_for_nu_over_full_range$nu_max, '           # *Maximum nu\n',
    curve_properties_for_nu_over_full_range$k_numin, '           # *Density at minimum nu\n',
    curve_properties_for_nu_over_full_range$nu_min, '           # *Minimum nu\n',
    curve_properties_for_nu_over_full_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_nu_over_full_range$n_troughs, '           # *No. of troughs\n',
    curve_properties_for_nu_over_full_range$nu_kjam, '           # Nu at jam density\n',
    curve_properties_for_nu_over_full_range$dnudk_kjam, '           # Gradient of the nu (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
cat('######################################################################################################################\n',
    '# ESTIMATED PROPERTIES OF THE FITTED "TAU" MODEL OVER THE DENSITY RANGE FROM ZERO TO "upper_density"\n',
    '# * = ESTIMATED FOR THE RESTRICTED DENSITY RANGE FROM ZERO TO THE ESTIMATED JAM DENSITY (IF DEFINED)\n',
    '######################################################################################################################\n',
    curve_properties_for_tau_over_full_range$tau_0, '           # Tau at zero density\n',
    curve_properties_for_tau_over_full_range$dtaudk_0, '           # Gradient of the tau (with respect to density) at zero density\n',
    curve_properties_for_tau_over_full_range$k_taumax, '           # *Density at maximum tau\n',
    curve_properties_for_tau_over_full_range$tau_max, '           # *Maximum tau\n',
    curve_properties_for_tau_over_full_range$k_taumin, '           # *Density at minimum tau\n',
    curve_properties_for_tau_over_full_range$tau_min, '           # *Minimum tau\n',
    curve_properties_for_tau_over_full_range$n_peaks, '           # *No. of peaks\n',
    curve_properties_for_tau_over_full_range$n_troughs, '           # *No. of troughs\n',
    curve_properties_for_tau_over_full_range$tau_kjam, '           # Tau at jam density\n',
    curve_properties_for_tau_over_full_range$dtaudk_kjam, '           # Gradient of the tau (with respect to density) at jam density\n',
    file = output_file, sep = '', append = TRUE)
}


################################################################################################################################################
plotA = function(traffic_data, reconstructed_model_fit, density_hi, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Fitted.Mu.For.XXXX.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create the plot object
ylo = min(0.0, traffic_data$V3, reconstructed_model_fit$mu)
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(limits = c(0.0, density_hi), expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(ylo, NA), expand = expand_scale(mult = 0.03)) +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_point(mapping = aes(x = V2, y = V3), data = traffic_data, colour = 'red', shape = 'circle small', size = 0.1) +
           geom_line(mapping = aes(x = V2, y = mu), data = reconstructed_model_fit, size = 0.5)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotB = function(traffic_data, density_hi, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Residuals.From.Mu.For.XXXX.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create the plot object
range_resid = range(traffic_data$V3 - traffic_data$fitted_values_mu)
ylo = min(0.0, range_resid[1])
yhi = max(0.0, range_resid[2])
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(limits = c(0.0, density_hi), expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(ylo, yhi), expand = expand_scale(mult = 0.03)) +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_point(mapping = aes(x = V2, y = V3 - fitted_values_mu), data = traffic_data, colour = 'red', shape = 'circle small', size = 0.1)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotC = function(traffic_data, reconstructed_model_fit, density_hi, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Percentiles.And.Mu.For.XXXX.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create the plot object
ylo = min(0.0, traffic_data$V3, reconstructed_model_fit$mu, reconstructed_model_fit$percentile_m3sig)
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(limits = c(0.0, density_hi), expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(ylo, NA), expand = expand_scale(mult = 0.03)) +
           geom_ribbon(mapping = aes(x = V2, ymin = percentile_m3sig, ymax = percentile_m2sig), data = reconstructed_model_fit, fill = 'grey90') +
           geom_ribbon(mapping = aes(x = V2, ymin = percentile_m2sig, ymax = percentile_m1sig), data = reconstructed_model_fit, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = V2, ymin = percentile_m1sig, ymax = percentile_p1sig), data = reconstructed_model_fit, fill = 'grey70') +
           geom_ribbon(mapping = aes(x = V2, ymin = percentile_p1sig, ymax = percentile_p2sig), data = reconstructed_model_fit, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = V2, ymin = percentile_p2sig, ymax = percentile_p3sig), data = reconstructed_model_fit, fill = 'grey90') +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_point(mapping = aes(x = V2, y = V3), data = traffic_data, colour = 'red', shape = 'circle small', size = 0.1) +
           geom_line(mapping = aes(x = V2, y = percentile_0sig), data = reconstructed_model_fit, colour = 'blue', linetype = 'dashed', size = 0.5) +
           geom_line(mapping = aes(x = V2, y = mu), data = reconstructed_model_fit, size = 0.5)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotD = function(traffic_data, density_hi, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Normalised.Quantile.Residuals.For.XXXX.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create the plot object
yhi = max(3.0, abs(traffic_data$normalised_quantile_residuals))
background_area = data.table(x = seq(from = 0.0, to = density_hi, length.out = 2))
background_area[, ylo := rep_len(-1.0, 2)]
background_area[, yhi := rep_len(1.0, 2)]
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(limits = c(0.0, density_hi), expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(-yhi, yhi), expand = expand_scale(mult = 0.03)) +
           geom_ribbon(mapping = aes(x = x, ymin = ylo - 2.0, ymax = yhi - 3.0), data = background_area, fill = 'grey90') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo - 1.0, ymax = yhi - 2.0), data = background_area, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo, ymax = yhi), data = background_area, fill = 'grey70') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo + 2.0, ymax = yhi + 1.0), data = background_area, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo + 3.0, ymax = yhi + 2.0), data = background_area, fill = 'grey90') +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_point(mapping = aes(x = V2, y = normalised_quantile_residuals), data = traffic_data, colour = 'red', shape = 'circle small', size = 0.1)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotE = function(traffic_data, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Normalised.Quantile.Residuals.Versus.Mu.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create the plot object
range_mu = range(traffic_data$fitted_values_mu)
xlo = min(0.0, range_mu[1])
xhi = max(0.0, range_mu[2])
yhi = max(3.0, abs(traffic_data$normalised_quantile_residuals))
background_area = data.table(x = seq(from = xlo, to = xhi, length.out = 2))
background_area[, ylo := rep_len(-1.0, 2)]
background_area[, yhi := rep_len(1.0, 2)]
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(limits = c(xlo, xhi), expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(-yhi, yhi), expand = expand_scale(mult = 0.03)) +
           geom_ribbon(mapping = aes(x = x, ymin = ylo - 2.0, ymax = yhi - 3.0), data = background_area, fill = 'grey90') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo - 1.0, ymax = yhi - 2.0), data = background_area, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo, ymax = yhi), data = background_area, fill = 'grey70') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo + 2.0, ymax = yhi + 1.0), data = background_area, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo + 3.0, ymax = yhi + 2.0), data = background_area, fill = 'grey90') +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_point(mapping = aes(x = fitted_values_mu, y = normalised_quantile_residuals), data = traffic_data, colour = 'red', shape = 'circle small', size = 0.1)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotF = function(traffic_data, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Normalised.Quantile.Residuals.Versus.Time.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Create the plot object
range_time = range(traffic_data$V1)
xlo = range_time[1]
xhi = range_time[2]
yhi = max(3.0, abs(traffic_data$normalised_quantile_residuals))
background_area = data.table(x = seq(from = xlo, to = xhi, length.out = 2))
background_area[, ylo := rep_len(-1.0, 2)]
background_area[, yhi := rep_len(1.0, 2)]
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(limits = c(xlo, xhi), expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(-yhi, yhi), expand = expand_scale(mult = 0.03)) +
           geom_ribbon(mapping = aes(x = x, ymin = ylo - 2.0, ymax = yhi - 3.0), data = background_area, fill = 'grey90') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo - 1.0, ymax = yhi - 2.0), data = background_area, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo, ymax = yhi), data = background_area, fill = 'grey70') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo + 2.0, ymax = yhi + 1.0), data = background_area, fill = 'grey80') +
           geom_ribbon(mapping = aes(x = x, ymin = ylo + 3.0, ymax = yhi + 2.0), data = background_area, fill = 'grey90') +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_point(mapping = aes(x = V1, y = normalised_quantile_residuals), data = traffic_data, colour = 'red', shape = 'circle small', size = 0.1)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotG = function(traffic_data, ntraffic_data, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Detrended.Normal.QQ.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>" (see
#              "FitFun.R" for details). The paper describing the detrended Normal quantile-quantile plot (or worm plot) is van Buuren & Fredriks
#              (2001, Statistics In Medicine, 20, 1259).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Prepare the data for plotting
zvals = qnorm(seq(from = 0.5/ntraffic_data, to = (ntraffic_data - 0.5)/ntraffic_data, length.out = ntraffic_data))
data_plot = data.table(zvals = zvals, detrended_nqr = sort(traffic_data$normalised_quantile_residuals) - zvals)

# Compute the 95% confidence interval
level = 0.95
fac = qnorm(0.5*(1.0 - level))
zseq = seq(from = min(-3.0, zvals[1]), to = max(3.0, zvals[ntraffic_data]), length.out = 1000)
pseq = pnorm(zseq)
ciseq = (fac/dnorm(zseq))*sqrt((pseq*(1 - pseq))/ntraffic_data)
ci_plot = data.table(zseq = zseq, ciseq = ciseq)

# Create the plot object
yhi = max(12.0/sqrt(ntraffic_data), abs(data_plot$detrended_nqr))
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str) +
           scale_x_continuous(expand = expand_scale(mult = 0.02)) +
           scale_y_continuous(limits = c(-yhi, yhi), expand = expand_scale(mult = 0.03)) +
           geom_hline(yintercept = 0, linetype = 'dotted') +
           geom_vline(xintercept = 0, linetype = 'dotted') +
           geom_line(mapping = aes(x = zseq, y = ciseq), data = ci_plot, linetype = 'dashed', size = 0.5) +
           geom_line(mapping = aes(x = zseq, y = -ciseq), data = ci_plot, linetype = 'dashed', size = 0.5) +
           geom_point(mapping = aes(x = zvals, y = detrended_nqr), data = data_plot, colour = 'red', shape = 'circle small', size = 0.1)

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
plotH = function(traffic_data, ntraffic_data, title_str, xlab_str, ylab_str, output_file) {

# Description: This function creates the plot "Plot.Of.Slotted.ACF.For.Normalised.Quantile.Residuals.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
#              (see "FitFun.R" for details).
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Determine the time lag bin size and compute the slotted auto-correlation function (slotted ACF). Note that for a small enough time lag bin size,
# the expected value of the slotted ACF in the zero time lag bin is unity since the slotted ACF is computed for a time-series of normalised quantile
# residuals (i.e. E(X^2) = 1 for X ~ N(0,1)).
tvec = traffic_data$V1
yvec = traffic_data$normalised_quantile_residuals
time_lag_bin_size = median(tvec[2:ntraffic_data] -  tvec[1:(ntraffic_data - 1)])
time_lag_nbins = 51
slotted_acf = compute_slotted_acf(tvec, yvec, time_lag_bin_size, time_lag_nbins)

# Create the plot object
plot_obj = ggplot() +
           theme_pubr(base_size = 16, border = TRUE) +
           theme(plot.title = element_text(size = 16, hjust = 0.5)) +
           ggtitle(title_str) +
           xlab(xlab_str) +
           ylab(ylab_str)

# Determine which values in the ACF data table should be plotted
selection_good = slotted_acf$acf_err > 0.0

# If there is at least one data point to be plotted
if (sum(selection_good) > 0) {

  # Update the plot object
  xlo = slotted_acf$tlag_bin_lo[1]
  xhi = slotted_acf$tlag_bin_hi[time_lag_nbins]
  slotted_acf = slotted_acf[selection_good]
  ylo = min(-1.0, slotted_acf$acf - slotted_acf$acf_err)
  yhi = max(1.0, slotted_acf$acf + slotted_acf$acf_err)
  plot_obj = plot_obj + scale_x_continuous(limits = c(xlo, xhi), expand = expand_scale(mult = 0.02)) +
                        scale_y_continuous(limits = c(ylo, yhi), expand = expand_scale(mult = 0.03)) +
                        geom_col(mapping = aes(x = tlag_bin_mid, y = acf), data = slotted_acf, width = time_lag_bin_size, fill = 'grey80') +
                        geom_hline(yintercept = 0, linetype = 'dotted') +
                        geom_errorbar(mapping = aes(x = tlag_bin_mid, ymin = acf - acf_err, ymax = acf + acf_err), data = slotted_acf, width = 0.3*time_lag_bin_size)
}

# Save the plot to the file "output_file"
ggsave(output_file, plot = plot_obj, scale = 2, width = 6.0, height = 4.0, units = 'in')
}


################################################################################################################################################
create_all_plots = function(traffic_data, ntraffic_data, data_max_density, upper_density, reconstructed_model_fit_selection, reconstructed_model_fit, ngrid,
                            fd_type, functional_form_model, noise_model, output_files) {

# Description: This function creates all of the plots for a GAMLSS model fit.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Define a useful variable
if (fd_type == 'Flow.Density') {
  fd_str = 'Flow'
} else if (fd_type == 'Speed.Density') {
  fd_str = 'Speed'
}

# Create a filtered version of the traffic data that does not include any normalised quantile residuals that are "-Inf" or "Inf"
selection = is.finite(traffic_data$normalised_quantile_residuals)
ntraffic_data_filtered = sum(selection)
if (ntraffic_data_filtered > 0) { traffic_data_filtered = traffic_data[selection] }

# Create the plot "Plot.Of.Fitted.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[4], '\n')
npts = nrow(reconstructed_model_fit_selection)
if (npts > 4000) {
  ind = ceiling((4000.0/npts)*seq(from = 1, to = npts))
  selection = rep_len(TRUE, npts)
  selection[2:(npts - 1)] = ind[2:(npts - 1)] != ind[1:(npts - 2)]
  reconstructed_model_fit_selection = reconstructed_model_fit_selection[selection]
}
title_str = paste0(fd_str, ' vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Data Density Range')
plotA(traffic_data, reconstructed_model_fit_selection, data_max_density, title_str, 'Density', fd_str, output_files[4])

# Create the plot "Plot.Of.Residuals.From.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[5], '\n')
title_str = paste0(fd_str, ' Residuals From Fitted Mu Curve vs Density : ', functional_form_model, ' : ', noise_model, ' : Data Density Range')
plotB(traffic_data, data_max_density, title_str, 'Density', paste0(fd_str, ' Residuals'), output_files[5])

# Create the plot "Plot.Of.Percentiles.And.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[6], '\n')
title_str = paste0(fd_str, ' vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Percentile Regions : Data Density Range')
plotC(traffic_data, reconstructed_model_fit_selection, data_max_density, title_str, 'Density', fd_str, output_files[6])

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
if (ntraffic_data_filtered > 0) {
  cat('Creating the plot:', output_files[7], '\n')
  title_str = paste0('Normalised Quantile Residuals (', fd_str, ') vs Density : ', functional_form_model, ' : ', noise_model, ' : Data Density Range')
  plotD(traffic_data_filtered, data_max_density, title_str, 'Density', paste0('Normalised Quantile Residuals (', fd_str, ')'), output_files[7])
} else {
  cat('WARNING - Cannot create the plot: ', output_files[7], '\n')
}

# Create the plot "Plot.Of.Fitted.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[8], '\n')
if (ngrid > 4000) {
  ind = ceiling((4000.0/ngrid)*seq(from = 1, to = ngrid))
  selection = rep_len(TRUE, ngrid)
  selection[2:(ngrid - 1)] = ind[2:(ngrid - 1)] != ind[1:(ngrid - 2)]
  reconstructed_model_fit = reconstructed_model_fit[selection]
}
title_str = paste0(fd_str, ' vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Full Density Range')
plotA(traffic_data, reconstructed_model_fit, upper_density, title_str, 'Density', fd_str, output_files[8])

# Create the plot "Plot.Of.Residuals.From.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[9], '\n')
title_str = paste0(fd_str, ' Residuals From Fitted Mu Curve vs Density : ', functional_form_model, ' : ', noise_model, ' : Full Density Range')
plotB(traffic_data, upper_density, title_str, 'Density', paste0(fd_str, ' Residuals'), output_files[9])

# Create the plot "Plot.Of.Percentiles.And.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[10], '\n')
title_str = paste0(fd_str, ' vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Percentile Regions : Full Density Range')
plotC(traffic_data, reconstructed_model_fit, upper_density, title_str, 'Density', fd_str, output_files[10])

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
if (ntraffic_data_filtered > 0) {
  cat('Creating the plot:', output_files[11], '\n')
  title_str = paste0('Normalised Quantile Residuals (', fd_str, ') vs Density : ', functional_form_model, ' : ', noise_model, ' : Full Density Range')
  plotD(traffic_data_filtered, upper_density, title_str, 'Density', paste0('Normalised Quantile Residuals (', fd_str, ')'), output_files[11])
} else {
  cat('WARNING - Cannot create the plot: ', output_files[11], '\n')
}

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.Versus.Mu.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
if (ntraffic_data_filtered > 0) {
  cat('Creating the plot:', output_files[12], '\n')
  title_str = paste0('Normalised Quantile Residuals (', fd_str, ') vs Fitted Mu Values : ', functional_form_model, ' : ', noise_model)
  plotE(traffic_data_filtered, title_str, 'Fitted Mu', paste0('Normalised Quantile Residuals (', fd_str, ')'), output_files[12])
} else {
  cat('WARNING - Cannot create the plot: ', output_files[12], '\n')
}

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.Versus.Time.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
if (ntraffic_data_filtered > 0) {
  cat('Creating the plot:', output_files[13], '\n')
  title_str = paste0('Normalised Quantile Residuals (', fd_str, ') vs Time : ', functional_form_model, ' : ', noise_model)
  plotF(traffic_data_filtered, title_str, 'Time', paste0('Normalised Quantile Residuals (', fd_str, ')'), output_files[13])
} else {
  cat('WARNING - Cannot create the plot: ', output_files[13], '\n')
}

# Create the plot "Plot.Of.Detrended.Normal.QQ.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
if (ntraffic_data_filtered > 0) {
  cat('Creating the plot:', output_files[14], '\n')
  title_str = paste0('Detrended Normal Q-Q Plot : ', functional_form_model, ' : ', noise_model, ' : 95% Confidence Interval')
  plotG(traffic_data_filtered, ntraffic_data_filtered, title_str, 'Standard Normal Quantiles (Units Of Sigma)', 'Detrended NQR Quantiles (Units Of Sigma)', output_files[14])
} else {
  cat('WARNING - Cannot create the plot: ', output_files[14], '\n')
}

# Create the plot "Plot.Of.Slotted.ACF.For.Normalised.Quantile.Residuals.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
if (ntraffic_data_filtered > 1) {
  cat('Creating the plot:', output_files[15], '\n')
  title_str = paste0('Slotted Auto-Correlation Function For Normalised Quantile Residuals : ', functional_form_model, ' : ', noise_model)
  plotH(traffic_data_filtered, ntraffic_data_filtered, title_str, 'Time Lag', 'Auto-Correlation Function', output_files[15])
} else {
  cat('WARNING - Cannot create the plot: ', output_files[15], '\n')
}
}


################################################################################################################################################
remove_file_list = function(file_list) {

# Description: This function deletes the files listed in "file_list".
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Delete the files in the file list "file_list"
for (file in file_list) {
  if (file.exists(file)) {
    file.remove(file)
  }
}
}
