fit_speed_density_with_DC1995Akjf_GaussSigCon = function(traffic_data, ngrid, upper_density, output_files) {

# Description: This function fits a GAMLSS model to the speed-density values in "traffic_data", and it is designed to be called directly from the R
#              script "FitFun.R". The model component for the functional form of the speed-density relationship is the Del Castillo model A with
#              fixed jam density (DC1995Akjf). The model component for the noise in the speed-density relationship is defined as independent
#              observations that follow a Gaussian distribution with constant variance (GaussSigCon).
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
k_jam = 1.0            # Fixed jam density (must be positive, and greater than or equal to the maximum observed density in the data)
par1_min = 0.0001      # Minimum acceptable value for the free parameter equivalent to (v_bw*k_jam)/(m*v_ff) (must be positive)
par1_step = 0.0001     # Step size for the free parameter equivalent to (v_bw*k_jam)/(m*v_ff) (must be positive)
par2_init = 1.0        # Initial value for the free parameter m (must be positive and greater than "par2_min")
par2_min = 0.0001      # Minimum acceptable value for the free parameter m (must be positive)
par2_step = 0.0001     # Step size for the free parameter m (must be positive)


# Define some useful variables
functional_form_model = 'DC1995Akjf'
noise_model = 'GaussSigCon'

# Report on the GAMLSS model and the data
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('The following GAMLSS model will be fit to the speed-density data:\n')
cat('\n')
cat('Model component for the functional form:\n')
cat('  Del Castillo\n')
cat('  Model A\n')
cat('  Fixed jam density (DC1995Akjf)\n')
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

  # Fit a Greenshields model with fixed jam density to estimate an initial value for v_ff
  { if (k_jam < data_max_density) {
      cat('ERROR - The jam density (fixed) is less than the maximum observed density...\n')
      q(save = 'no', status = 1)
    }
    init_model_obj = gamlss(V3 ~ 0 + I(1.0 - (V2/k_jam)), sigma.formula = ~ 1, family = NO(), data = traffic_data)
    if (init_model_obj$converged != TRUE) {
      cat('ERROR - The initial fit of a Greenshields model with fixed jam density did not converge...\n')
      q(save = 'no', status = 1)
    }
    v_ff_init = init_model_obj$mu.coefficients[1]

    # Fit a Greenberg model with fixed jam density to estimate an initial value for (v_bw*k_jam)/(m*v_ff)
    init_model_obj = gamlss(V3 ~ 0 + log(k_jam/V2), sigma.formula = ~ 1, family = NO(), data = traffic_data)
    if (init_model_obj$converged != TRUE) {
      cat('ERROR - The initial fit of a Greenberg model with fixed jam density did not converge...\n')
      q(save = 'no', status = 1)
    }
    par1_init = max((init_model_obj$mu.coefficients[1]*k_jam)/(par2_init*v_ff_init), par1_min + par1_step)

    # Perform the intermediate fits
    k_jam_use = data.frame(k_jam_use = k_jam)
    model_formula = quote(gamlss(V3 ~ 0 + I(1.0 - exp(1.0 - ((1.0 + (p[1]*((1.0/V2) - (1.0/k_jam_use))))^p[2]))), sigma.formula = ~ 1, family = NO()))
    attach(k_jam_use)
    attach(traffic_data)
    optim_obj = try(find.hyper(model = model_formula, parameters = c(par1_init, par2_init), k = 0.0, steps = c(par1_step, par2_step), lower = c(par1_min, par2_min),
                               maxit = 500))
    if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
    if (optim_obj$convergence != 0) {
      par1_min_use = data.frame(par1_min_use = par1_min)
      par2_min_use = data.frame(par2_min_use = par2_min)
      model_formula = quote(gamlss(V3 ~ 0 + I(1.0 - exp(1.0 - ((1.0 + ((par1_min_use + abs(p[1]))*((1.0/V2) - (1.0/k_jam_use))))^(par2_min_use + abs(p[2]))))),
                                   sigma.formula = ~ 1, family = NO()))
      attach(par1_min_use)
      attach(par2_min_use)
      optim_obj = try(find.hyper(model = model_formula, parameters = c(par1_init - par1_min, par2_init - par2_min), k = 0.0, steps = c(par1_step, par2_step),
                                 method = 'Nelder-Mead', maxit = 500))
      if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
      detach(par2_min_use)
      detach(par1_min_use)
      detach(traffic_data)
      detach(k_jam_use)
      if (optim_obj$convergence != 0) {
        cat('ERROR - The intermediate fits did not converge...\n')
        q(save = 'no', status = 1)
      }
      par1 = par1_min + abs(optim_obj$par[1])
      par2 = par2_min + abs(optim_obj$par[2])
    } else {
      detach(traffic_data)
      detach(k_jam_use)
      par1 = optim_obj$par[1]
      par2 = optim_obj$par[2]
    }

    # Perform the final fit
    model_obj = gamlss(V3 ~ 0 + I(1.0 - exp(1.0 - ((1.0 + (par1*((1.0/V2) - (1.0/k_jam))))^par2))), sigma.formula = ~ 1, family = NO(), data = traffic_data)
    if (model_obj$converged != TRUE) {
      cat('ERROR - The final fit did not converge...\n')
      q(save = 'no', status = 1)
    }
    model_obj$mu.df = model_obj$mu.df + 2
    model_obj$df.fit = model_obj$df.fit + 2
    model_obj$df.residual = model_obj$df.residual - 2
    model_obj$aic = model_obj$aic + 4.0
    model_obj$sbc = model_obj$sbc + 2.0*log(ntraffic_data) },
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
  { reconstructed_model_fit = data.table(V2 = seq(from = 0.0, to = min(upper_density, k_jam), length.out = ngrid))
    predicted_values_for_mu = double(length = ngrid)
    predicted_values_for_mu[2:ngrid] = predict(model_obj, what = 'mu', newdata = reconstructed_model_fit[2:ngrid], type = 'response', data = traffic_data)
    predicted_values_for_mu[1] = model_obj$mu.coefficients[1]
    predicted_values_for_sigma = predict(model_obj, what = 'sigma', newdata = reconstructed_model_fit, type = 'response', data = traffic_data)
    if (!all(is.finite(predicted_values_for_mu))) {
      cat('ERROR - The reconstructed fitted model for "mu" includes at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (!all(is.finite(predicted_values_for_sigma))) {
      cat('ERROR - The reconstructed fitted model for "sigma" includes at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(predicted_values_for_sigma <= 0.0)) {
      cat('ERROR - The reconstructed fitted model for "sigma" includes at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    reconstructed_model_fit[, mu := predicted_values_for_mu]
    reconstructed_model_fit[, sigma := predicted_values_for_sigma]
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
    dvdk_0 = NA     # Can be computed analytically - not yet implemented
    k_crit = NA
    k_vmax = NA
    q_cap = NA
    v_max = NA
    v_bw = NA
    dvdk_kjam = NA
    if (model_obj$mu.coefficients[1] > 0.0) {
      v_ff = model_obj$mu.coefficients[1]
#      k_crit =     # Can be computed numerically - not yet implemented
      k_vmax = 0.0
#      q_cap =      # Can be computed once k_crit is available - not yet implemented
      v_max = v_ff
      v_bw = par1*par2*(v_ff/k_jam)
      dvdk_kjam = -v_bw/k_jam
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
cat('  v_ff:                 ', sprintf('%.8g', model_obj$mu.coefficients[1]), '\n')
cat('  (v_bw*k_jam)/(m*v_ff):', sprintf('%.8g', par1), '\n')
cat('  m:                    ', sprintf('%.8g', par2), '\n')
cat('  sigma_con:            ', sprintf('%.8g', exp(model_obj$sigma.coefficients[1])), '\n')

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
        sprintf('%.8g', par1), '           # (v_bw*k_jam)/(m*v_ff)\n',
        sprintf('%.8g', par2), '           # m\n',
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
