fit_flow_density_with_FF_GCV = function(data, ngrid, upper_density, output_files) {

# Description: This function fits a GAMLSS model to the flow-density values in "data", and it is designed to be called directly from the R script
#              "FitFun.R". The model component for the functional form of the flow-density relationship is the free-flow model (FF). The model
#              component for the noise in the flow-density relationship is defined as independent observations that follow a Gaussian distribution
#              with constant variance (GCV).
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
#   NONE


# Define some useful variables
functional_form_model = 'FF'
noise_model = 'GCV'

# Report on the GAMLSS model and the data
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('The following GAMLSS model will be fit to the flow-density data:\n')
cat('\n')
cat('Model component for the functional form:\n')
cat('  Free-flow\n')
cat('\n')
cat('Model component for the noise:\n')
cat('  Independent observations\n')
cat('  Gaussian distribution\n')
cat('  Constant variance\n')
cat('\n')
cat('Data properties:\n')
tryCatch(
  { ndata = nrow(data)
    data_range_density = range(data$V2)
    data_min_density = data_range_density[1]
    data_max_density = data_range_density[2]
    data_range_flow = range(data$V3)
    data_min_flow = data_range_flow[1]
    data_max_flow = data_range_flow[2] },
  error = function(cond) { cat('ERROR - Failed to determine the data properties...\n')
                           q(save = 'no', status = 1) }
)
cat('  No. of flow-density measurement pairs (Ndat):', ndata, '\n')
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
  { model_obj = gamlss(V3 ~ 0 + V2, sigma.formula = ~ 1, family = NO(), data = data) },
  error = function(cond) { cat('ERROR - Failed to fit the GAMLSS model...\n')
                           q(save = 'no', status = 1) }
)

# Check that the model fit converged
if (model_obj$converged != TRUE) {
  cat('ERROR - The fit did not converge...\n')
  q(save = 'no', status = 1)
}

# Store the predicted values for the model at the density values in the data, along with the normalised quantile residuals, in the data table
cat('\n')
cat('Storing the predicted values for the model, along with the normalised quantile residuals, in the data table...\n')
tryCatch(
  { if (!all(is.finite(model_obj$mu.fv))) {
      cat('ERROR - The predicted values for "mu" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (!all(is.finite(model_obj$sigma.fv))) {
      cat('ERROR - The predicted values for "sigma" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (!all(is.finite(model_obj$residuals))) {
      cat('ERROR - The normalised quantile residuals include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(model_obj$sigma.fv <= 0.0)) {
      cat('ERROR - The predicted values for "sigma" at the density values in the data include at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    data[, fitted_values_mu := model_obj$mu.fv]
    data[, fitted_values_sigma := model_obj$sigma.fv]
    data[, fitted_values_nu := double(length = ndata)]
    data[, fitted_values_tau := rep_len(3.0, ndata)]
    data[, normalised_quantile_residuals := model_obj$residuals] },
  error = function(cond) { cat('ERROR - Failed to store the predicted values for the model and the normalised quantile residuals...\n')
                           q(save = 'no', status = 1) }
)

# Reconstruct the fitted model over the density range from zero to "upper_density"
cat('Reconstructing the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { reconstructed_model_fit = data.table(V2 = seq(from = 0.0, to = upper_density, length.out = ngrid))
    predicted_values = predictAll(model_obj, newdata = reconstructed_model_fit, type = 'response', data = data)
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
    reconstructed_model_fit[, tau := rep_len(3.0, ngrid)] },
  error = function(cond) { cat('ERROR - Failed to reconstruct the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Construct percentile curves for the fitted model over the density range from zero to "upper_density"
cat('Constructing percentile curves for the fitted model over the density range from 0 to', upper_density, '...\n')
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
    v_ff = model_obj$mu.coefficients[1]
    dvdk_0 = 0.0
    k_crit = NA
    k_vmax = NA
    q_cap = NA
    v_max = NA
    k_jam = NA
    v_bw = NA
    dvdk_kjam = NA },
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
cat('  Global deviance (-2 ln L):   ', gdev, '\n')
cat('  AIC (-2 ln L + 2 Npar):      ', aic, '\n')
cat('  BIC (-2 ln L + Npar ln Ndat):', bic, '\n')
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
cat('Fitted model parameters (see the accompanying paper by Bramich, Menendez & Ambuhl for details):\n')
cat('  v_ff:     ', model_obj$mu.coefficients[1], '\n')
cat('  sigma_con:', exp(model_obj$sigma.coefficients[1]), '\n')

# Write out the fit summary file "Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt"
cat('\n')
cat('Writing out the fit summary file:    ', output_files[1], '\n')
tryCatch(
  { write_fit_summary(output_files[1], 'Flow.Density', ndata, data_min_density, data_max_density, data_min_flow, data_max_flow,
                      npar_mu, npar_sigma, npar_nu, npar_tau, npar_all, gdev, aic, bic,
                      q_0, v_ff, dvdk_0, k_crit, k_vmax, q_cap, v_max, k_jam, v_bw, dvdk_kjam,
                      curve_properties_for_mu_over_data_range, curve_properties_for_sigma_over_data_range,
                      curve_properties_for_nu_over_data_range, curve_properties_for_tau_over_data_range,
                      curve_properties_for_mu_over_full_range, curve_properties_for_sigma_over_full_range,
                      curve_properties_for_nu_over_full_range, curve_properties_for_tau_over_full_range)
    cat('######################################################################################################################\n',
        '# FITTED MODEL PARAMETERS (SEE THE ACCOMPANYING PAPER BY BRAMICH, MENENDEZ & AMBUHL FOR DETAILS)\n',
        '# N.B: FITTED COEFFICIENTS FOR ANY NON-PARAMETRIC SMOOTHING FUNCTIONS IN THE MODEL ARE NOT REPORTED HERE\n',
        '######################################################################################################################\n',
        model_obj$mu.coefficients[1], '           # v_ff\n',
        exp(model_obj$sigma.coefficients[1]), '           # sigma_con\n',
        file = output_files[1], sep = '', append = TRUE)
    cat('######################################################################################################################\n',
        '# FIT SUMMARY AS PROVIDED BY THE GAMLSS SOFTWARE\n',
        '######################################################################################################################\n',
        file = output_files[1], sep = '', append = TRUE)
    sink(file = output_files[1], append = TRUE)
    summary(model_obj)
    sink() },
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
        'Normalised Quantile Residual\n', file = output_files[3])
    write.table(data, file = output_files[3], append = TRUE, quote = FALSE, row.names = FALSE, col.names = FALSE) },
  error = function(cond) { cat('ERROR - Failed to write out the fit predictions file...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# If no plots should be produced, then return the GAMLSS model fit object
if (length(output_files) == 3) {
  cat('\n')
  cat('>-----------------------------------------------------------------------------<\n')
  return(model_obj)
}

# Create the plot "Plot.Of.Fitted.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('\n')
cat('Creating the plot:', output_files[4], '\n')
tryCatch(
  { npts = nrow(reconstructed_model_fit_selection)
    if (npts > 4000) {
      ind = ceiling((4000.0/npts)*seq(from = 1, to = npts))
      selection = rep_len(TRUE, npts)
      selection[2:(npts - 1)] = ind[2:(npts - 1)] != ind[1:(npts - 2)]
      reconstructed_model_fit_selection = reconstructed_model_fit_selection[selection]
    }
    title_str = paste0('Flow vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Data Density Range')
    plotA(data, reconstructed_model_fit_selection, title_str, 'Density', 'Flow', output_files[4]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Residuals.From.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[5], '\n')
tryCatch(
  { title_str = paste0('Flow Residuals From Fitted Mu Curve vs Density : ', functional_form_model, ' : ', noise_model, ' : Data Density Range')
    plotB(data, data_max_density, title_str, 'Density', 'Flow Residuals', output_files[5]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Percentiles.And.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[6], '\n')
tryCatch(
  { title_str = paste0('Flow vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Percentile Regions : Data Density Range')
    plotC(data, reconstructed_model_fit_selection, title_str, 'Density', 'Flow', output_files[6]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[7], '\n')
tryCatch(
  { title_str = paste0('Normalised Quantile Residuals (Flow) vs Density : ', functional_form_model, ' : ', noise_model, ' : Data Density Range')
    plotD(data, data_max_density, title_str, 'Density', 'Normalised Quantile Residuals (Flow)', output_files[7]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Fitted.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[8], '\n')
tryCatch(
  { if (ngrid > 4000) {
      ind = ceiling((4000.0/ngrid)*seq(from = 1, to = ngrid))
      selection = rep_len(TRUE, ngrid)
      selection[2:(ngrid - 1)] = ind[2:(ngrid - 1)] != ind[1:(ngrid - 2)]
      reconstructed_model_fit = reconstructed_model_fit[selection]
    }
    title_str = paste0('Flow vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Full Density Range')
    plotA(data, reconstructed_model_fit, title_str, 'Density', 'Flow', output_files[8]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Residuals.From.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[9], '\n')
tryCatch(
  { title_str = paste0('Flow Residuals From Fitted Mu Curve vs Density : ', functional_form_model, ' : ', noise_model, ' : Full Density Range')
    plotB(data, upper_density, title_str, 'Density', 'Flow Residuals', output_files[9]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Percentiles.And.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[10], '\n')
tryCatch(
  { title_str = paste0('Flow vs Density : ', functional_form_model, ' : ', noise_model, ' : Fitted Mu Curve : Percentile Regions : Full Density Range')
    plotC(data, reconstructed_model_fit, title_str, 'Density', 'Flow', output_files[10]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[11], '\n')
tryCatch(
  { title_str = paste0('Normalised Quantile Residuals (Flow) vs Density : ', functional_form_model, ' : ', noise_model, ' : Full Density Range')
    plotD(data, upper_density, title_str, 'Density', 'Normalised Quantile Residuals (Flow)', output_files[11]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.Versus.Mu.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[12], '\n')
tryCatch(
  { title_str = paste0('Normalised Quantile Residuals (Flow) vs Fitted Mu Values : ', functional_form_model, ' : ', noise_model)
    plotE(data, title_str, 'Fitted Mu', 'Normalised Quantile Residuals (Flow)', output_files[12]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Normalised.Quantile.Residuals.Versus.Time.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[13], '\n')
tryCatch(
  { title_str = paste0('Normalised Quantile Residuals (Flow) vs Time : ', functional_form_model, ' : ', noise_model)
    plotF(data, title_str, 'Time', 'Normalised Quantile Residuals (Flow)', output_files[13]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Detrended.Normal.QQ.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[14], '\n')
tryCatch(
  { title_str = paste0('Detrended Normal Q-Q Plot : ', functional_form_model, ' : ', noise_model, ' : 95% Confidence Interval')
    plotG(data, ndata, title_str, 'Theoretical Quantiles (Units Of Sigma)', 'Deviation From Theoretical Quantiles (Units Of Sigma)', output_files[14]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Create the plot "Plot.Of.Slotted.ACF.For.Normalised.Quantile.Residuals.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>"
cat('Creating the plot:', output_files[15], '\n')
tryCatch(
  { title_str = paste0('Slotted Auto-Correlation Function For Normalised Quantile Residuals : ', functional_form_model, ' : ', noise_model)
    plotH(data, ndata, title_str, 'Time Lag', 'Auto-Correlation Function', output_files[15]) },
  error = function(cond) { cat('ERROR - Failed to create the plot...\n')
                           remove_file_list(output_files)
                           q(save = 'no', status = 1) }
)

# Return the GAMLSS model fit object
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
return(model_obj)
}
