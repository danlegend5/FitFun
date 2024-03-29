fit_speed_density_with_DK1966A_GaussSigCon = function(traffic_data, ngrid, upper_density, output_files) {

# Description: This function fits a GAMLSS model to the speed-density values in "traffic_data", and it is designed to be called directly from the R
#              script "FitFun.R". The model component for the functional form of the speed-density relationship is the Drake multi-regime model A
#              (DK1966A). The model component for the noise in the speed-density relationship is defined as independent observations that follow a
#              Gaussian distribution with constant variance (GaussSigCon).
#                The input parameters "ngrid" and "upper_density" are used to define an equally spaced grid of "ngrid" density values ranging from
#              zero to "upper_density". The function employs this density grid to reconstruct the fitted model at the grid points for use in plots
#              and for estimating certain properties of the fitted model that are not directly accessible from the fitted parameter values.
#                The function creates various output files "output_files" including diagnostic plots (see "FitFun.R" for details). If the function
#              finishes successfully, then it returns the corresponding GAMLSS model fit object.
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)
#
# Configuration Parameters:
#
par1_nsteps = 50     # Number of steps to be used for the initial profiling of the free parameter k_b (must be greater than or equal to 2)
par1_nrefine = 7     # Number of refinement iterations to be performed for fitting the free parameter k_b (must be greater than or equal to 1)


# Define some useful variables
functional_form_model = 'DK1966A'
noise_model = 'GaussSigCon'

# Report on the GAMLSS model and the data
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('The following GAMLSS model will be fit to the speed-density data:\n')
cat('\n')
cat('Model component for the functional form:\n')
cat('  Drake\n')
cat('  Multi-regime\n')
cat('  Model A (DK1966A)\n')
cat('\n')
cat('Model component for the noise:\n')
cat('  Independent observations\n')
cat('  Gaussian distribution\n')
cat('  Constant variance (GaussSigCon)\n')
cat('\n')
cat('Data properties:\n')
tryCatch(
  { ntraffic_data = nrow(traffic_data)
    data_range_density = range(traffic_data$V2)
    data_min_density = data_range_density[1]
    data_max_density = data_range_density[2]
    data_range_speed = range(traffic_data$V3)
    data_min_speed = data_range_speed[1]
    data_max_speed = data_range_speed[2] },
  error = function(cond) { cat('ERROR - Failed to determine the data properties...\n')
                           q(save = 'no', status = 1) }
)
cat('  No. of speed-density measurement pairs (Ndat):', ntraffic_data, '\n')
cat('  Minimum density in the data:                  ', sprintf('%.8g', data_min_density), '\n')
cat('  Maximum density in the data:                  ', sprintf('%.8g', data_max_density), '\n')
cat('  Minimum speed in the data:                    ', sprintf('%.8g', data_min_speed), '\n')
cat('  Maximum speed in the data:                    ', sprintf('%.8g', data_max_speed), '\n')
cat('\n')
cat('Model reconstruction:\n')
tryCatch(
  { grid_density_step = upper_density/(ngrid - 1) },
  error = function(cond) { cat('ERROR - Failed to compute the grid density step...\n')
                           q(save = 'no', status = 1) }
)
cat('  No. of density grid points:', ngrid, '\n')
cat('  Grid lower density:         0\n')
cat('  Grid upper density:        ', sprintf('%.8g', upper_density), '\n')
cat('  Grid density step:         ', sprintf('%.8g', grid_density_step), '\n')

# Fit the GAMLSS model to the data
cat('\n')
cat('Fitting the GAMLSS model...\n')
tryCatch(

  # Perform the fits for the initial profiling of k_b
  { k_b_vec = sort(traffic_data$V2)
    k_b_lo = k_b_vec[2]
    k_b_hi = k_b_vec[ntraffic_data - 1]
    k_b_vec = seq(from = k_b_lo, to = k_b_hi, length.out = par1_nsteps + 2)
    k_b_vec = k_b_vec[2:(par1_nsteps + 1)]
    gdev_vec = double(length = par1_nsteps)
    for (i in 1:par1_nsteps) {
      curr_k_b = k_b_vec[i]
      model_obj = gamlss(V3 ~ 0 + ifelse(V2 > curr_k_b, 0.0, 1.0) + ifelse(V2 > curr_k_b, 0.0, V2) + ifelse(V2 > curr_k_b, 1.0, 0.0) + ifelse(V2 > curr_k_b, V2, 0.0),
                         sigma.formula = ~ 1, family = NO(), data = traffic_data)
      if (model_obj$converged != TRUE) {
        cat('ERROR - The initial profiling fits did not converge...\n')
        q(save = 'no', status = 1)
      }
      gdev_vec[i] = model_obj$G.deviance
    }
    min_gdev = min(gdev_vec)
    ind_min_gdev = which(gdev_vec == min_gdev)
    ind_min_gdev = ind_min_gdev[ceiling(0.5*length(ind_min_gdev))]

    # Iteratively refine the best fitting model
    for (j in 1:par1_nrefine) {
      ind_lo = max(1, ind_min_gdev - 3)
      ind_hi = min(par1_nsteps, ind_min_gdev + 3)
      curr_par1_nsteps = 2*(ind_hi - ind_lo) + 1
      curr_k_b_vec = double(length = curr_par1_nsteps)
      curr_gdev_vec = rep_len(NA, curr_par1_nsteps)
      for (i in ind_lo:ind_hi) {
        tmp_ind = 2*(i - ind_lo) + 1
        curr_k_b_vec[tmp_ind] = k_b_vec[i]
        curr_gdev_vec[tmp_ind] = gdev_vec[i]
      }
      for (i in 1:(ind_hi - ind_lo)) {
        tmp_ind = 2*i
        curr_k_b_vec[tmp_ind] = 0.5*(curr_k_b_vec[tmp_ind - 1] + curr_k_b_vec[tmp_ind + 1])
      }
      for (i in 1:curr_par1_nsteps) {
        if (!is.na(curr_gdev_vec[i])) { next }
        curr_k_b = curr_k_b_vec[i]
        model_obj = gamlss(V3 ~ 0 + ifelse(V2 > curr_k_b, 0.0, 1.0) + ifelse(V2 > curr_k_b, 0.0, V2) + ifelse(V2 > curr_k_b, 1.0, 0.0) + ifelse(V2 > curr_k_b, V2, 0.0),
                           sigma.formula = ~ 1, family = NO(), data = traffic_data)
        if (model_obj$converged != TRUE) {
          cat('ERROR - The refining fits did not converge...\n')
          q(save = 'no', status = 1)
        }
        curr_gdev_vec[i] = model_obj$G.deviance
      }
      curr_min_gdev = min(curr_gdev_vec)
      curr_ind_min_gdev = which(curr_gdev_vec == curr_min_gdev)
      curr_ind_min_gdev = curr_ind_min_gdev[ceiling(0.5*length(curr_ind_min_gdev))]
      par1_nsteps = curr_par1_nsteps
      k_b_vec = curr_k_b_vec
      gdev_vec = curr_gdev_vec
      ind_min_gdev = curr_ind_min_gdev
    }
    par1 = k_b_vec[ind_min_gdev]

    # Perform the final fit
    model_obj = gamlss(V3 ~ 0 + ifelse(V2 > par1, 0.0, 1.0) + ifelse(V2 > par1, 0.0, V2) + ifelse(V2 > par1, 1.0, 0.0) + ifelse(V2 > par1, V2, 0.0),
                       sigma.formula = ~ 1, family = NO(), data = traffic_data)
    if (model_obj$converged != TRUE) {
      cat('ERROR - The final fit did not converge...\n')
      q(save = 'no', status = 1)
    }
    model_obj$mu.df = model_obj$mu.df + 1
    model_obj$df.fit = model_obj$df.fit + 1
    model_obj$df.residual = model_obj$df.residual - 1
    model_obj$aic = model_obj$aic + 2.0
    model_obj$sbc = model_obj$sbc + log(ntraffic_data) },
  error = function(cond) { cat('ERROR - Failed to fit the GAMLSS model...\n')
                           q(save = 'no', status = 1) }
)

# Store the predicted values for the fitted model at the density values in the data in the data table
cat('\n')
cat('Storing the predicted values for the fitted model in the data table...\n')
tryCatch(
  { if (!all(is.finite(model_obj$mu.fv))) {
      cat('ERROR - The predicted values for "mu" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (!all(is.finite(model_obj$sigma.fv))) {
      cat('ERROR - The predicted values for "sigma" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(model_obj$sigma.fv <= 0.0)) {
      cat('ERROR - The predicted values for "sigma" at the density values in the data include at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    traffic_data[, fitted_values_mu := model_obj$mu.fv]
    traffic_data[, fitted_values_sigma := model_obj$sigma.fv]
    traffic_data[, fitted_values_nu := double(length = ntraffic_data)]
    traffic_data[, fitted_values_tau := double(length = ntraffic_data)] },
  error = function(cond) { cat('ERROR - Failed to store the predicted values for the fitted model...\n')
                           q(save = 'no', status = 1) }
)

# Compute the normalised quantile residuals and store them in the data table. Note that the normalised quantile residuals may include some "-Inf"
# or "Inf" values for particularly bad outliers.
cat('Computing the normalised quantile residuals...\n')
tryCatch(
  { cumulative_probs_lower = pNO(traffic_data$V3, mu = model_obj$mu.fv, sigma = model_obj$sigma.fv)
    cumulative_probs_upper = pNO(traffic_data$V3, mu = model_obj$mu.fv, sigma = model_obj$sigma.fv, lower.tail = FALSE)
    traffic_data[, normalised_quantile_residuals := calculate_normalised_quantile_residuals(cumulative_probs_lower, cumulative_probs_upper)] },
  error = function(cond) { cat('ERROR - Failed to compute the normalised quantile residuals...\n')
                           q(save = 'no', status = 1) }
)

# Compute the percentiles for the fitted model at the density values in the data and store them in the data table
cat('Computing the percentiles for the fitted model at the density values in the data...\n')
tryCatch(
  { traffic_data[, percentile_m3sig := qNO(pNO(-3.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)]
    traffic_data[, percentile_m2sig := qNO(pNO(-2.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)]
    traffic_data[, percentile_m1sig := qNO(pNO(-1.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)]
    traffic_data[, percentile_0sig := qNO(0.5, mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)]
    traffic_data[, percentile_p1sig := qNO(pNO(1.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)]
    traffic_data[, percentile_p2sig := qNO(pNO(2.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)]
    traffic_data[, percentile_p3sig := qNO(pNO(3.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma)] },
  error = function(cond) { cat('ERROR - Failed to compute the percentiles for the fitted model at the density values in the data...\n')
                           q(save = 'no', status = 1) }
)

# Compute a set of distributional measures for the fitted model at the density values in the data and store them in the data table
cat('Computing a set of distributional measures for the fitted model at the density values in the data...\n')
tryCatch(
  { traffic_data[, mean := traffic_data$fitted_values_mu]
    traffic_data[, median := traffic_data$fitted_values_mu]
    traffic_data[, mode := traffic_data$fitted_values_mu]
    traffic_data[, standard_deviation := traffic_data$fitted_values_sigma]
    traffic_data[, moment_skewness := traffic_data$fitted_values_nu]
    traffic_data[, moment_excess_kurtosis := traffic_data$fitted_values_tau] },
  error = function(cond) { cat('ERROR - Failed to compute a set of distributional measures for the fitted model at the density values in the data...\n')
                           q(save = 'no', status = 1) }
)

# Reconstruct the fitted model over the density range from zero to "upper_density"
cat('Reconstructing the fitted model over the density range from 0 to', sprintf('%.8g', upper_density), '...\n')
tryCatch(
  { reconstructed_model_fit = data.table(V2 = seq(from = 0.0, to = upper_density, length.out = ngrid))
    predicted_values = predictAll(model_obj, newdata = reconstructed_model_fit, type = 'response', data = traffic_data)
    if (!all(is.finite(predicted_values$mu))) {
      cat('ERROR - The reconstructed fitted model for "mu" includes at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (!all(is.finite(predicted_values$sigma))) {
      cat('ERROR - The reconstructed fitted model for "sigma" includes at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(predicted_values$sigma <= 0.0)) {
      cat('ERROR - The reconstructed fitted model for "sigma" includes at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    reconstructed_model_fit[, mu := predicted_values$mu]
    reconstructed_model_fit[, sigma := predicted_values$sigma]
    reconstructed_model_fit[, nu := double(length = ngrid)]
    reconstructed_model_fit[, tau := double(length = ngrid)] },
  error = function(cond) { cat('ERROR - Failed to reconstruct the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Construct percentile curves for the fitted model over the density range from zero to "upper_density"
cat('Constructing percentile curves for the fitted model over the density range from 0 to', sprintf('%.8g', upper_density), '...\n')
tryCatch(
  { reconstructed_model_fit[, percentile_m3sig := qNO(pNO(-3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, percentile_m2sig := qNO(pNO(-2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, percentile_m1sig := qNO(pNO(-1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, percentile_0sig := qNO(0.5, mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, percentile_p1sig := qNO(pNO(1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, percentile_p2sig := qNO(pNO(2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, percentile_p3sig := qNO(pNO(3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)] },
  error = function(cond) { cat('ERROR - Failed to construct percentile curves for the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Construct curves of a set of distributional measures for the fitted model over the density range from zero to "upper_density"
cat('Constructing curves of a set of distributional measures for the fitted model over the density range from 0 to', sprintf('%.8g', upper_density), '...\n')
tryCatch(
  { reconstructed_model_fit[, mean := reconstructed_model_fit$mu]
    reconstructed_model_fit[, median := reconstructed_model_fit$mu]
    reconstructed_model_fit[, mode := reconstructed_model_fit$mu]
    reconstructed_model_fit[, standard_deviation := reconstructed_model_fit$sigma]
    reconstructed_model_fit[, moment_skewness := reconstructed_model_fit$nu]
    reconstructed_model_fit[, moment_excess_kurtosis := reconstructed_model_fit$tau] },
  error = function(cond) { cat('ERROR - Failed to construct curves of a set of distributional measures for the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Estimate useful properties of the fitted model over the density range from zero to the maximum observed density using the reconstruction
cat('Estimating useful properties of the fitted model using the reconstruction...\n')
tryCatch(
  { selection = reconstructed_model_fit$V2 < (data_max_density + grid_density_step)
    reconstructed_model_fit_selection = reconstructed_model_fit[selection]
    curve_properties_for_mu_over_data_range = get_curve_properties_for_mu(reconstructed_model_fit_selection, 'Speed.Density')
    curve_properties_for_sigma_over_data_range = get_curve_properties_for_sigma(reconstructed_model_fit_selection, curve_properties_for_mu_over_data_range)
    curve_properties_for_nu_over_data_range = get_curve_properties_for_nu(reconstructed_model_fit_selection, curve_properties_for_mu_over_data_range)
    curve_properties_for_tau_over_data_range = get_curve_properties_for_tau(reconstructed_model_fit_selection, curve_properties_for_mu_over_data_range) },
  error = function(cond) { cat('ERROR - Failed to estimate useful properties of the fitted model over the density range from zero to the maximum observed density...\n')
                           q(save = 'no', status = 1) }
)

# Estimate useful properties of the fitted model over the density range from zero to "upper_density" using the reconstruction
tryCatch(
  { curve_properties_for_mu_over_full_range = get_curve_properties_for_mu(reconstructed_model_fit, 'Speed.Density')
    curve_properties_for_sigma_over_full_range = get_curve_properties_for_sigma(reconstructed_model_fit, curve_properties_for_mu_over_full_range)
    curve_properties_for_nu_over_full_range = get_curve_properties_for_nu(reconstructed_model_fit, curve_properties_for_mu_over_full_range)
    curve_properties_for_tau_over_full_range = get_curve_properties_for_tau(reconstructed_model_fit, curve_properties_for_mu_over_full_range) },
  error = function(cond) { cat('ERROR - Failed to estimate useful properties of the fitted model over the density range from zero to "upper_density"...\n')
                           q(save = 'no', status = 1) }
)

# Extract information from the model fit object for the fit summary
cat('Extracting information from the model fit object for the fit summary...\n')
tryCatch(
  { npar_mu = model_obj$mu.df
    npar_sigma = model_obj$sigma.df
    npar_nu = 0
    npar_tau = 0
    npar_all = model_obj$df.fit
    gdev = model_obj$G.deviance
    aic = model_obj$aic
    bic = model_obj$sbc },
  error = function(cond) { cat('ERROR - Failed to extract information from the model fit object for the fit summary...\n')
                           q(save = 'no', status = 1) }
)

# Where possible, extract physical parameter values from the model fit object for the fit summary
tryCatch(
  { q_0 = 0.0
    v_ff = NA
    dvdk_0 = model_obj$mu.coefficients[2]
    k_crit = NA     # Can be computed analytically - not yet implemented
    k_vmax = NA     # Can be computed analytically - not yet implemented
    q_cap = NA      # Can be computed once k_crit is available - not yet implemented
    v_max = NA      # Can be computed once k_vmax is available - not yet implemented
    k_jam = NA
    v_bw = NA
    dvdk_kjam = NA
    if (model_obj$mu.coefficients[1] > 0.0) { v_ff = model_obj$mu.coefficients[1] }
    if (model_obj$mu.coefficients[3] > 0.0) {
      if (model_obj$mu.coefficients[4] < 0.0) {
        if (-model_obj$mu.coefficients[3]/model_obj$mu.coefficients[4] > par1) {
          v_bw = model_obj$mu.coefficients[3]
          dvdk_kjam = model_obj$mu.coefficients[4]
          k_jam = -v_bw/dvdk_kjam
        }
      }
    } },
  error = function(cond) { cat('ERROR - Failed to extract physical parameter values from the model fit object for the fit summary...\n')
                           q(save = 'no', status = 1) }
)

# Report a basic fit summary
cat('\n')
cat('Fit summary (abridged):\n')
cat('\n')
cat('Model parameter counts:\n')
cat('  No. of free parameters (mu):        ', sprintf('%.6f', npar_mu), '\n')
cat('  No. of free parameters (sigma):     ', sprintf('%.6f', npar_sigma), '\n')
cat('  No. of free parameters (nu):        ', sprintf('%.6f', npar_nu), '\n')
cat('  No. of free parameters (tau):       ', sprintf('%.6f', npar_tau), '\n')
cat('  Total no. of free parameters (Npar):', sprintf('%.6f', npar_all), '\n')
cat('\n')
cat('Fit quality:\n')
cat('  Global deviance [ -2*ln(Lmax) ]:    ', sprintf('%.4f', gdev), '\n')
cat('  AIC [ -2*ln(Lmax) + 2*Npar ]:       ', sprintf('%.4f', aic), '\n')
cat('  BIC [ -2*ln(Lmax) + Npar*ln(Ndat) ]:', sprintf('%.4f', bic), '\n')
cat('\n')
cat('Fitted physical parameters (where available):\n')
if (!is.na(q_0)) { cat('  Flow at zero density:                                  ', sprintf('%.8g', q_0), '\n') }
if (!is.na(v_ff)) { cat('  Free-flow speed:                                       ', sprintf('%.8g', v_ff), '\n') }
if (!is.na(dvdk_0)) { cat('  Gradient of the speed (w.r.t. density) at zero density:', sprintf('%.8g', dvdk_0), '\n') }
if (!is.na(k_crit)) { cat('  Critical density:                                      ', sprintf('%.8g', k_crit), '\n') }
if (!is.na(k_vmax)) { cat('  Density at maximum speed:                              ', sprintf('%.8g', k_vmax), '\n') }
if (!is.na(q_cap)) { cat('  Capacity:                                              ', sprintf('%.8g', q_cap), '\n') }
if (!is.na(v_max)) { cat('  Maximum speed:                                         ', sprintf('%.8g', v_max), '\n') }
if (!is.na(k_jam)) { cat('  Jam density:                                           ', sprintf('%.8g', k_jam), '\n') }
if (!is.na(v_bw)) { cat('  Back-propagating wave speed at jam density:            ', sprintf('%.8g', v_bw), '\n') }
if (!is.na(dvdk_kjam)) { cat('  Gradient of the speed (w.r.t. density) at jam density: ', sprintf('%.8g', dvdk_kjam), '\n') }
cat('\n')
cat('Fitted model parameters (see the accompanying papers by Bramich, Menendez & Ambuhl for details):\n')
cat('  v_ff:      ', sprintf('%.8g', model_obj$mu.coefficients[1]), '\n')
cat('  c:         ', sprintf('%.8g', -model_obj$mu.coefficients[2]), '\n')
cat('  k_b:       ', sprintf('%.8g', par1), '\n')
cat('  v_bw:      ', sprintf('%.8g', model_obj$mu.coefficients[3]), '\n')
cat('  v_bw/k_jam:', sprintf('%.8g', -model_obj$mu.coefficients[4]), '\n')
cat('  sigma_con: ', sprintf('%.8g', exp(model_obj$sigma.coefficients[1])), '\n')

# Write out the fit summary file "Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt"
cat('\n')
cat('Writing out the fit summary file:    ', output_files[1], '\n')
tryCatch(
  { write_fit_summary(output_files[1], 'Speed.Density', ntraffic_data, data_min_density, data_max_density, data_min_speed, data_max_speed,
                      npar_mu, npar_sigma, npar_nu, npar_tau, npar_all, gdev, aic, bic,
                      q_0, v_ff, dvdk_0, k_crit, k_vmax, q_cap, v_max, k_jam, v_bw, dvdk_kjam,
                      curve_properties_for_mu_over_data_range, curve_properties_for_sigma_over_data_range,
                      curve_properties_for_nu_over_data_range, curve_properties_for_tau_over_data_range,
                      curve_properties_for_mu_over_full_range, curve_properties_for_sigma_over_full_range,
                      curve_properties_for_nu_over_full_range, curve_properties_for_tau_over_full_range)
    cat('######################################################################################################################\n',
        '# FITTED MODEL PARAMETERS (SEE THE ACCOMPANYING PAPERS BY BRAMICH, MENENDEZ & AMBUHL FOR DETAILS)\n',
        '# N.B: FITTED COEFFICIENTS FOR ANY NON-PARAMETRIC SMOOTHING FUNCTIONS IN THE MODEL ARE NOT REPORTED HERE\n',
        '######################################################################################################################\n',
        sprintf('%.8g', model_obj$mu.coefficients[1]), '           # v_ff\n',
        sprintf('%.8g', -model_obj$mu.coefficients[2]), '           # c\n',
        sprintf('%.8g', par1), '           # k_b\n',
        sprintf('%.8g', model_obj$mu.coefficients[3]), '           # v_bw\n',
        sprintf('%.8g', -model_obj$mu.coefficients[4]), '           # v_bw/k_jam\n',
        sprintf('%.8g', exp(model_obj$sigma.coefficients[1])), '           # sigma_con\n',
        file = output_files[1], sep = '', append = TRUE) },
  error = function(cond) { cat('ERROR - Failed to write out the fit summary file...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Write out the fit curves file "Fit.Curves.<fd_type>.<functional_form_model>.<noise_model>.txt"
cat('Writing out the fit curves file:     ', output_files[2], '\n')
tryCatch(
  { cat('# Density : Mu : Sigma : Nu : Tau : 0.135 Percentile (Corresponding To -3*Sigma In A Normal Distribution) : 2.28 Percentile (Corresponding To -2*Sigma In A',
        'Normal Distribution) : 15.87 Percentile (Corresponding To -1*Sigma In A Normal Distribution) : 50.00 Percentile (Median) : 84.13 Percentile (Corresponding',
        'To 1*Sigma In A Normal Distribution) : 97.72 Percentile (Corresponding To 2*Sigma In A Normal Distribution) : 99.865 Percentile (Corresponding To 3*Sigma',
        'In A Normal Distribution) : Mean : Median : Mode : Standard Deviation : Moment Skewness : Moment Excess Kurtosis\n', file = output_files[2])
    write.table(reconstructed_model_fit, file = output_files[2], append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) },
  error = function(cond) { cat('ERROR - Failed to write out the fit curves file...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Write out the fit predictions file "Fit.Predictions.<fd_type>.<functional_form_model>.<noise_model>.txt"
cat('Writing out the fit predictions file:', output_files[3], '\n')
tryCatch(
  { cat('# Data Column 1 : Data Column 2 : Data Column 3 : Fitted Value For Mu : Fitted Value For Sigma : Fitted Value For Nu : Fitted Value For Tau :',
        'Normalised Quantile Residual ("-Inf" Or "Inf" Values May Be Present) : 0.135 Percentile (Corresponding To -3*Sigma In A Normal Distribution) :',
        '2.28 Percentile (Corresponding To -2*Sigma In A Normal Distribution) : 15.87 Percentile (Corresponding To -1*Sigma In A Normal Distribution) :',
        '50.00 Percentile (Median) : 84.13 Percentile (Corresponding To 1*Sigma In A Normal Distribution) : 97.72 Percentile (Corresponding To 2*Sigma',
        'In A Normal Distribution) : 99.865 Percentile (Corresponding To 3*Sigma In A Normal Distribution) : Mean : Median : Mode : Standard Deviation :',
        'Moment Skewness : Moment Excess Kurtosis\n', file = output_files[3])
    write.table(traffic_data, file = output_files[3], append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) },
  error = function(cond) { cat('ERROR - Failed to write out the fit predictions file...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# If required, then create the plots for the GAMLSS model fit
if (length(output_files) > 3) {
  cat('\n')
  cat('Creating the plots for the GAMLSS model fit...\n')
  tryCatch(
    { create_all_plots(traffic_data, ntraffic_data, data_max_density, upper_density, reconstructed_model_fit_selection, reconstructed_model_fit, ngrid,
                       'Speed.Density', functional_form_model, noise_model, output_files) },
    error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                             remove_file_list(output_files)
                             q(save = 'no', status = 1) }
  )
}

# Return the GAMLSS model fit object
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
return(model_obj)
}
