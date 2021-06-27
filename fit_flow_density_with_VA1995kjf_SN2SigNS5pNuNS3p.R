fit_flow_density_with_VA1995kjf_SN2SigNS5pNuNS3p = function(traffic_data, ngrid, upper_density, output_files) {

# Description: This function fits a GAMLSS model to the flow-density values in "traffic_data", and it is designed to be called directly from the R
#              script "FitFun.R". The model component for the functional form of the flow-density relationship is the Van Aerde model with fixed
#              jam density (VA1995kjf). The model component for the noise in the flow-density relationship is defined as independent observations
#              that follow a Skew Normal Type II distribution. The density dependence of the log of the scale parameter and the log of the skewness
#              parameter is modelled using natural cubic splines with five and three effective free parameters, respectively (SN2SigNS5pNuNS3p).
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
k_jam = 1.0             # Fixed jam density (must be positive)
psi_min = 0.0001        # Minimum acceptable value for the free parameter psi (must be positive)
psi_step = 0.0001       # Step size for the free parameter psi (must be positive)
omega_min = 0.0001      # Minimum acceptable value for the free parameter omega (must be positive)
omega_step = 0.0001     # Step size for the free parameter omega (must be positive)
inner_ccrit = 0.05      # Convergence criterion for the inner iteration of the GAMLSS fitting algorithm
inner_ncyc = 10         # Maximum number of cycles of the inner iteration of the GAMLSS fitting algorithm
outer_ccrit = 0.05      # Convergence criterion for the outer iteration of the GAMLSS fitting algorithm
outer_ncyc = 1000       # Maximum number of cycles of the outer iteration of the GAMLSS fitting algorithm


# Define some useful variables
functional_form_model = 'VA1995kjf'
noise_model = 'SN2SigNS5pNuNS3p'

# Report on the GAMLSS model and the data
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('The following GAMLSS model will be fit to the flow-density data:\n')
cat('\n')
cat('Model component for the functional form:\n')
cat('  Van Aerde\n')
cat('  Fixed jam density (VA1995kjf)\n')
cat('\n')
cat('Model component for the noise:\n')
cat('  Independent observations\n')
cat('  Skew Normal Type II distribution\n')
cat('  Scale is a smooth function of density\n')
cat('  Skewness is a smooth function of density (SN2SigNS5pNuNS3p)\n')
cat('\n')
cat('Data properties:\n')
tryCatch(
  { ntraffic_data = nrow(traffic_data)
    data_range_density = range(traffic_data$V2)
    data_min_density = data_range_density[1]
    data_max_density = data_range_density[2]
    data_range_flow = range(traffic_data$V3)
    data_min_flow = data_range_flow[1]
    data_max_flow = data_range_flow[2] },
  error = function(cond) { cat('ERROR - Failed to determine the data properties...\n')
                           q(save = 'no', status = 1) }
)
cat('  No. of flow-density measurement pairs (Ndat):', ntraffic_data, '\n')
cat('  Minimum density in the data:                 ', data_min_density, '\n')
cat('  Maximum density in the data:                 ', data_max_density, '\n')
cat('  Minimum flow in the data:                    ', data_min_flow, '\n')
cat('  Maximum flow in the data:                    ', data_max_flow, '\n')
cat('\n')
cat('Model reconstruction:\n')
tryCatch(
  { grid_density_step = upper_density/(ngrid - 1) },
  error = function(cond) { cat('ERROR - Failed to compute the grid density step...\n')
                           q(save = 'no', status = 1) }
)
cat('  No. of density grid points:', ngrid, '\n')
cat('  Grid lower density:         0\n')
cat('  Grid upper density:        ', upper_density, '\n')
cat('  Grid density step:         ', grid_density_step, '\n')

# Fit the GAMLSS model to the data
cat('\n')
cat('Fitting the GAMLSS model...\n')
tryCatch(

  # Fit a Greenshields model with fixed jam density to estimate an initial value for v_ff
  { init_model_obj = gamlss(V3 ~ 0 + I(V2*(1.0 - (V2/k_jam))), sigma.formula = ~ 1, family = NO(), data = traffic_data)
    if (init_model_obj$converged != TRUE) {
      cat('ERROR - The initial fit of a Greenshields model with fixed jam density did not converge...\n')
      q(save = 'no', status = 1)
    }
    v_ff_init = init_model_obj$mu.coefficients[1]

    # Fit a Greenberg model with fixed jam density to estimate an initial value for v_bw
    init_model_obj = gamlss(V3 ~ 0 + I(V2*log(k_jam/V2)), sigma.formula = ~ 1, family = NO(), data = traffic_data)
    if (init_model_obj$converged != TRUE) {
      cat('ERROR - The initial fit of a Greenberg model with fixed jam density did not converge...\n')
      q(save = 'no', status = 1)
    }
    v_bw_init = init_model_obj$mu.coefficients[1]

    # Compute an initial value for the Van Aerde parameter c3
    c3_init = max(((exp(1.0)/v_bw_init) - (4.0/v_ff_init))/k_jam, 0.0)

    # Compute initial values for the Bramich parameters psi and omega
    psi_init = max(1.0/k_jam, psi_min + psi_step)
    omega_init = max(c3_init*v_ff_init, omega_min + omega_step)

    # Perform the intermediate fits
    inv_k_jam = data.frame(inv_k_jam = 1.0/k_jam)
    inner_ccrit_use = data.frame(inner_ccrit_use = inner_ccrit)
    inner_ncyc_use = data.frame(inner_ncyc_use = inner_ncyc)
    outer_ccrit_use = data.frame(outer_ccrit_use = outer_ccrit)
    outer_ncyc_use = data.frame(outer_ncyc_use = outer_ncyc)
    model_formula = quote(gamlss(V3 ~ 0 + I(1.0 - ((inv_k_jam - p[1] - p[2])*V2) - sqrt(((((inv_k_jam - p[1] + p[2])*V2) - 1.0)^2) + (4.0*p[1]*p[2]*(V2^2)))),
                                 sigma.formula = ~ ns(V2, df = 4), nu.formula = ~ ns(V2, df = 2), family = SN2(), c.crit = outer_ccrit_use, n.cyc = outer_ncyc_use,
                                 i.control = glim.control(cc = inner_ccrit_use, cyc = inner_ncyc_use)))
    attach(inv_k_jam)
    attach(traffic_data)
    attach(inner_ccrit_use)
    attach(inner_ncyc_use)
    attach(outer_ccrit_use)
    attach(outer_ncyc_use)
    optim_obj = try(find.hyper(model = model_formula, parameters = c(psi_init, omega_init), k = 0.0, steps = c(psi_step, omega_step), lower = c(psi_min, omega_min),
                               maxit = 500))
    if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
    if (optim_obj$convergence != 0) {
      psi_min_use = data.frame(psi_min_use = psi_min)
      omega_min_use = data.frame(omega_min_use = omega_min)
      model_formula = quote(gamlss(V3 ~ 0 + I(1.0 - ((inv_k_jam - (psi_min_use + abs(p[1])) - (omega_min_use + abs(p[2])))*V2) - sqrt(((((inv_k_jam - (psi_min_use + abs(p[1])) + (omega_min_use + abs(p[2])))*V2) - 1.0)^2) + (4.0*(psi_min_use + abs(p[1]))*(omega_min_use + abs(p[2]))*(V2^2)))),
                                   sigma.formula = ~ ns(V2, df = 4), nu.formula = ~ ns(V2, df = 2), family = SN2(), c.crit = outer_ccrit_use, n.cyc = outer_ncyc_use,
                                   i.control = glim.control(cc = inner_ccrit_use, cyc = inner_ncyc_use)))
      attach(psi_min_use)
      attach(omega_min_use)
      optim_obj = try(find.hyper(model = model_formula, parameters = c(psi_init - psi_min, omega_init - omega_min), k = 0.0, steps = c(psi_step, omega_step),
                                 method = 'Nelder-Mead', maxit = 500))
      if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
      detach(omega_min_use)
      detach(psi_min_use)
      detach(outer_ncyc_use)
      detach(outer_ccrit_use)
      detach(inner_ncyc_use)
      detach(inner_ccrit_use)
      detach(traffic_data)
      detach(inv_k_jam)
      if (optim_obj$convergence != 0) {
        cat('ERROR - The intermediate fits did not converge...\n')
        q(save = 'no', status = 1)
      }
      psi = psi_min + abs(optim_obj$par[1])
      omega = omega_min + abs(optim_obj$par[2])
    } else {
      detach(outer_ncyc_use)
      detach(outer_ccrit_use)
      detach(inner_ncyc_use)
      detach(inner_ccrit_use)
      detach(traffic_data)
      detach(inv_k_jam)
      psi = optim_obj$par[1]
      omega = optim_obj$par[2]
    }

    # Perform the final fit
    model_obj = gamlss(V3 ~ 0 + I(1.0 - (((1.0/k_jam) - psi - omega)*V2) - sqrt((((((1.0/k_jam) - psi + omega)*V2) - 1.0)^2) + (4.0*psi*omega*(V2^2)))),
                       sigma.formula = ~ ns(V2, df = 4), nu.formula = ~ ns(V2, df = 2), family = SN2(), data = traffic_data,
                       c.crit = outer_ccrit, n.cyc = outer_ncyc, i.control = glim.control(cc = inner_ccrit, cyc = inner_ncyc))
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

# Store the predicted values for the model at the density values in the data in the data table
cat('\n')
cat('Storing the predicted values for the model in the data table...\n')
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
    if (!all(is.finite(model_obj$nu.fv))) {
      cat('ERROR - The predicted values for "nu" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(model_obj$nu.fv <= 0.0)) {
      cat('ERROR - The predicted values for "nu" at the density values in the data include at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    traffic_data[, fitted_values_mu := model_obj$mu.fv]
    traffic_data[, fitted_values_sigma := model_obj$sigma.fv]
    traffic_data[, fitted_values_nu := model_obj$nu.fv]
    traffic_data[, fitted_values_tau := double(length = ntraffic_data)] },
  error = function(cond) { cat('ERROR - Failed to store the predicted values for the model...\n')
                           q(save = 'no', status = 1) }
)

# Compute the normalised quantile residuals and store them in the data table. Note that the normalised quantile residuals may include some "-Inf"
# or "Inf" values for particularly bad outliers.
cat('Computing the normalised quantile residuals and storing them in the data table...\n')
tryCatch(
  { cumulative_probs_lower = pSN2(traffic_data$V3, mu = model_obj$mu.fv, sigma = model_obj$sigma.fv, nu = model_obj$nu.fv)
    cumulative_probs_upper = pSN2(traffic_data$V3, mu = model_obj$mu.fv, sigma = model_obj$sigma.fv, nu = model_obj$nu.fv, lower.tail = FALSE)
    traffic_data[, normalised_quantile_residuals := calculate_normalised_quantile_residuals(cumulative_probs_lower, cumulative_probs_upper)] },
  error = function(cond) { cat('ERROR - Failed to compute and store the normalised quantile residuals...\n')
                           q(save = 'no', status = 1) }
)

# Reconstruct the fitted model over the density range from zero to "upper_density"
cat('Reconstructing the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { reconstructed_model_fit = data.table(V2 = seq(from = 0.0, to = upper_density, length.out = ngrid))
    predicted_values = predictAll(model_obj, newdata = reconstructed_model_fit, type = 'response', data = traffic_data)
    predicted_values$sigma = fix_out_of_data_curves(reconstructed_model_fit$V2, predicted_values$sigma, data_min_density, data_max_density)
    predicted_values$nu = fix_out_of_data_curves(reconstructed_model_fit$V2, predicted_values$nu, data_min_density, data_max_density)
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
    if (!all(is.finite(predicted_values$nu))) {
      cat('ERROR - The reconstructed fitted model for "nu" includes at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(predicted_values$nu <= 0.0)) {
      cat('ERROR - The reconstructed fitted model for "nu" includes at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    reconstructed_model_fit[, mu := predicted_values$mu]
    reconstructed_model_fit[, sigma := predicted_values$sigma]
    reconstructed_model_fit[, nu := predicted_values$nu]
    reconstructed_model_fit[, tau := double(length = ngrid)] },
  error = function(cond) { cat('ERROR - Failed to reconstruct the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Construct percentile curves for the fitted model over the density range from zero to "upper_density"
cat('Constructing percentile curves for the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { reconstructed_model_fit[, percentile_m3sig := qSN2(pNO(-3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)]
    reconstructed_model_fit[, percentile_m2sig := qSN2(pNO(-2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)]
    reconstructed_model_fit[, percentile_m1sig := qSN2(pNO(-1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)]
    reconstructed_model_fit[, percentile_0sig := qSN2(0.5, mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)]
    reconstructed_model_fit[, percentile_p1sig := qSN2(pNO(1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)]
    reconstructed_model_fit[, percentile_p2sig := qSN2(pNO(2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)]
    reconstructed_model_fit[, percentile_p3sig := qSN2(pNO(3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma, nu = reconstructed_model_fit$nu)] },
  error = function(cond) { cat('ERROR - Failed to construct percentile curves for the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Estimate useful properties of the fitted model over the density range from zero to the maximum observed density using the reconstruction
cat('Estimating useful properties of the fitted model using the reconstruction...\n')
tryCatch(
  { selection = reconstructed_model_fit$V2 < (data_max_density + grid_density_step)
    reconstructed_model_fit_selection = reconstructed_model_fit[selection]
    curve_properties_for_mu_over_data_range = get_curve_properties_for_mu(reconstructed_model_fit_selection, 'Flow.Density')
    curve_properties_for_sigma_over_data_range = get_curve_properties_for_sigma(reconstructed_model_fit_selection, curve_properties_for_mu_over_data_range)
    curve_properties_for_nu_over_data_range = get_curve_properties_for_nu(reconstructed_model_fit_selection, curve_properties_for_mu_over_data_range)
    curve_properties_for_tau_over_data_range = get_curve_properties_for_tau(reconstructed_model_fit_selection, curve_properties_for_mu_over_data_range) },
  error = function(cond) { cat('ERROR - Failed to estimate useful properties of the fitted model over the density range from zero to the maximum observed density...\n')
                           q(save = 'no', status = 1) }
)

# Estimate useful properties of the fitted model over the density range from zero to "upper_density" using the reconstruction
tryCatch(
  { curve_properties_for_mu_over_full_range = get_curve_properties_for_mu(reconstructed_model_fit, 'Flow.Density')
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
    npar_nu = model_obj$nu.df
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
    k_crit = NA     # Can be computed analytically - not yet implemented
    k_vmax = NA     # Can be computed analytically - not yet implemented
    q_cap = NA      # Can be computed once k_crit is available - not yet implemented
    v_max = NA      # Can be computed once k_vmax is available - not yet implemented
    v_bw = NA
    dvdk_kjam = NA
    if (model_obj$mu.coefficients[1] > 0.0) {
      if (omega > 0.0) {
        v_ff = 2.0*model_obj$mu.coefficients[1]*omega
        v_bw = v_ff/(k_jam*(psi + omega))
        dvdk_kjam = -v_bw/k_jam
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
cat('  No. of free parameters (mu):        ', npar_mu, '\n')
cat('  No. of free parameters (sigma):     ', npar_sigma, '\n')
cat('  No. of free parameters (nu):        ', npar_nu, '\n')
cat('  No. of free parameters (tau):       ', npar_tau, '\n')
cat('  Total no. of free parameters (Npar):', npar_all, '\n')
cat('\n')
cat('Fit quality:\n')
cat('  Global deviance [ -2*ln(Lmax) ]:    ', format(gdev, nsmall = 4), '\n')
cat('  AIC [ -2*ln(Lmax) + 2*Npar ]:       ', format(aic, nsmall = 4), '\n')
cat('  BIC [ -2*ln(Lmax) + Npar*ln(Ndat) ]:', format(bic, nsmall = 4), '\n')
cat('\n')
cat('Fitted physical parameters (where available):\n')
if (!is.na(q_0)) { cat('  Flow at zero density:                                  ', q_0, '\n') }
if (!is.na(v_ff)) { cat('  Free-flow speed:                                       ', v_ff, '\n') }
if (!is.na(dvdk_0)) { cat('  Gradient of the speed (w.r.t. density) at zero density:', dvdk_0, '\n') }
if (!is.na(k_crit)) { cat('  Critical density:                                      ', k_crit, '\n') }
if (!is.na(k_vmax)) { cat('  Density at maximum speed:                              ', k_vmax, '\n') }
if (!is.na(q_cap)) { cat('  Capacity:                                              ', q_cap, '\n') }
if (!is.na(v_max)) { cat('  Maximum speed:                                         ', v_max, '\n') }
if (!is.na(k_jam)) { cat('  Jam density:                                           ', k_jam, '\n') }
if (!is.na(v_bw)) { cat('  Back-propagating wave speed at jam density:            ', v_bw, '\n') }
if (!is.na(dvdk_kjam)) { cat('  Gradient of the speed (w.r.t. density) at jam density: ', dvdk_kjam, '\n') }
cat('\n')
cat('Fitted model parameters (see the accompanying papers by Bramich, Menendez & Ambuhl for details):\n')
cat('  alpha:', model_obj$mu.coefficients[1], '\n')
cat('  psi:  ', psi, '\n')
cat('  omega:', omega, '\n')

# Write out the fit summary file "Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt"
cat('\n')
cat('Writing out the fit summary file:    ', output_files[1], '\n')
tryCatch(
  { write_fit_summary(output_files[1], 'Flow.Density', ntraffic_data, data_min_density, data_max_density, data_min_flow, data_max_flow,
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
        model_obj$mu.coefficients[1], '           # alpha\n',
        psi, '           # psi\n',
        omega, '           # omega\n',
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
        'In A Normal Distribution)\n', file = output_files[2])
    write.table(reconstructed_model_fit, file = output_files[2], append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) },
  error = function(cond) { cat('ERROR - Failed to write out the fit curves file...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Write out the fit predictions file "Fit.Predictions.<fd_type>.<functional_form_model>.<noise_model>.txt"
cat('Writing out the fit predictions file:', output_files[3], '\n')
tryCatch(
  { cat('# Data Column 1 : Data Column 2 : Data Column 3 : Fitted Value For Mu : Fitted Value For Sigma : Fitted Value For Nu : Fitted Value For Tau :',
        'Normalised Quantile Residual ("-Inf" Or "Inf" Values May Be Present)\n', file = output_files[3])
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
                       'Flow.Density', functional_form_model, noise_model, output_files) },
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
