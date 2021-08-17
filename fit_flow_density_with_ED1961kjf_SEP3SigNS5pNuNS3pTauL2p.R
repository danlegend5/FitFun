fit_flow_density_with_ED1961kjf_SEP3SigNS5pNuNS3pTauL2p = function(traffic_data, ngrid, upper_density, output_files) {

# Description: This function fits a GAMLSS model to the flow-density values in "traffic_data", and it is designed to be called directly from the R
#              script "FitFun.R". The model component for the functional form of the flow-density relationship is the Edie multi-regime model with
#              fixed jam density (ED1961kjf). The model component for the noise in the flow-density relationship is defined as independent
#              observations that follow a Skew Exponential Power Type III distribution. The density dependence of the log of the scale parameter
#              and the log of the skewness parameter is modelled using natural cubic splines with five and three effective free parameters,
#              respectively. The density dependence of the log of the kurtosis parameter is modelled using a straight line function (i.e. an
#              intercept and a gradient as the two free parameters; SEP3SigNS5pNuNS3pTauL2p).
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
par1_step = 0.0001     # Step size for the free parameter equivalent to -1/k_crit (must be positive)
par2_nsteps = 50       # Number of steps to be used for the initial profiling of the free parameter k_b (must be greater than or equal to 2)
par2_nrefine = 7       # Number of refinement iterations to be performed for fitting the free parameter k_b (must be greater than or equal to 1)
inner_ccrit = 0.05     # Convergence criterion for the inner iteration of the GAMLSS fitting algorithm
inner_ncyc = 10        # Maximum number of cycles of the inner iteration of the GAMLSS fitting algorithm
outer_ccrit = 0.05     # Convergence criterion for the outer iteration of the GAMLSS fitting algorithm
outer_ncyc = 1000      # Maximum number of cycles of the outer iteration of the GAMLSS fitting algorithm


# Define some useful variables
functional_form_model = 'ED1961kjf'
noise_model = 'SEP3SigNS5pNuNS3pTauL2p'

# Report on the GAMLSS model and the data
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('The following GAMLSS model will be fit to the flow-density data:\n')
cat('\n')
cat('Model component for the functional form:\n')
cat('  Edie\n')
cat('  Multi-regime\n')
cat('  Fixed jam density (ED1961kjf)\n')
cat('\n')
cat('Model component for the noise:\n')
cat('  Independent observations\n')
cat('  Skew Exponential Power Type III distribution\n')
cat('  Scale is a smooth function of density\n')
cat('  Skewness is a smooth function of density\n')
cat('  Kurtosis is a smooth function of density (SEP3SigNS5pNuNS3pTauL2p)\n')
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

  # Fit an Underwood model A to estimate an initial value for -1/k_crit
  { if (k_jam < data_max_density) {
      cat('ERROR - The jam density (fixed) is less than the maximum observed density...\n')
      q(save = 'no', status = 1)
    }
    init_model_obj = gamlss(V3 ~ 1 + offset(log(V2)) + V2, sigma.formula = ~ 1, family = NO(mu.link = 'log'), data = traffic_data)
    if (init_model_obj$converged != TRUE) {
      cat('ERROR - The initial fit of an Underwood model A did not converge...\n')
      q(save = 'no', status = 1)
    }
    par1_init = init_model_obj$mu.coefficients[2]

    # Perform the fits for the initial profiling of k_b
    k_b_vec = sort(traffic_data$V2)
    k_b_lo = k_b_vec[2]
    k_b_hi = k_b_vec[ntraffic_data - 1]
    k_b_vec = seq(from = k_b_lo, to = k_b_hi, length.out = par2_nsteps + 2)
    k_b_vec = k_b_vec[2:(par2_nsteps + 1)]
    par1_vec = double(length = par2_nsteps)
    gdev_vec = double(length = par2_nsteps)
    k_jam_use = data.frame(k_jam_use = k_jam)
    inner_ccrit_use = data.frame(inner_ccrit_use = inner_ccrit)
    inner_ncyc_use = data.frame(inner_ncyc_use = inner_ncyc)
    outer_ccrit_use = data.frame(outer_ccrit_use = outer_ccrit)
    outer_ncyc_use = data.frame(outer_ncyc_use = outer_ncyc)
    model_formula = quote(gamlss(V3 ~ 0 + ifelse(V2 > curr_k_b, 0.0, V2*exp(p[1]*V2)) + ifelse(V2 > curr_k_b, V2*log(k_jam_use/V2), 0.0),
                                 sigma.formula = ~ ns(V2, df = 4), nu.formula = ~ ns(V2, df = 2), tau.formula = ~ 1 + V2, family = SEP3(),
                                 c.crit = outer_ccrit_use, n.cyc = outer_ncyc_use, i.control = glim.control(cc = inner_ccrit_use, cyc = inner_ncyc_use)))
    for (i in 1:par2_nsteps) {
      curr_k_b = data.frame(curr_k_b = k_b_vec[i])
      attach(k_jam_use)
      attach(curr_k_b)
      attach(traffic_data)
      attach(inner_ccrit_use)
      attach(inner_ncyc_use)
      attach(outer_ccrit_use)
      attach(outer_ncyc_use)
      optim_obj = try(find.hyper(model = model_formula, parameters = c(par1_init), k = 0.0, steps = c(par1_step), maxit = 500))
      if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
      if (optim_obj$convergence != 0) {
        par1_range = max(10.0*abs(par1_init), 1000.0*par1_step)
        par1_min = par1_init - par1_range
        par1_max = par1_init + par1_range
        optim_obj = try(find.hyper(model = model_formula, parameters = c(par1_init), k = 0.0, steps = c(par1_step), lower = c(par1_min), upper = c(par1_max),
                                   method = 'Brent', maxit = 500))
        if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
        if (optim_obj$convergence != 0) {
          cat('ERROR - The initial profiling fits did not converge...\n')
          detach(outer_ncyc_use)
          detach(outer_ccrit_use)
          detach(inner_ncyc_use)
          detach(inner_ccrit_use)
          detach(traffic_data)
          detach(curr_k_b)
          detach(k_jam_use)
          q(save = 'no', status = 1)
        }
        if ((optim_obj$par[1] <= (par1_min + par1_step)) || (optim_obj$par[1] >= (par1_max - par1_step))) {
          cat('ERROR - The initial profiling fits did not converge (parameter limit reached)...\n')
          detach(outer_ncyc_use)
          detach(outer_ccrit_use)
          detach(inner_ncyc_use)
          detach(inner_ccrit_use)
          detach(traffic_data)
          detach(curr_k_b)
          detach(k_jam_use)
          q(save = 'no', status = 1)
        }
      }
      detach(outer_ncyc_use)
      detach(outer_ccrit_use)
      detach(inner_ncyc_use)
      detach(inner_ccrit_use)
      detach(traffic_data)
      detach(curr_k_b)
      detach(k_jam_use)
      par1_vec[i] = optim_obj$par[1]
      gdev_vec[i] = optim_obj$value
      par1_init = par1_vec[i]
    }
    min_gdev = min(gdev_vec)
    ind_min_gdev = which(gdev_vec == min_gdev)
    ind_min_gdev = ind_min_gdev[ceiling(0.5*length(ind_min_gdev))]

    # Iteratively refine the best fitting model
    for (j in 1:par2_nrefine) {
      ind_lo = max(1, ind_min_gdev - 3)
      ind_hi = min(par2_nsteps, ind_min_gdev + 3)
      curr_par2_nsteps = 2*(ind_hi - ind_lo) + 1
      curr_k_b_vec = double(length = curr_par2_nsteps)
      curr_par1_vec = double(length = curr_par2_nsteps)
      curr_gdev_vec = rep_len(NA, curr_par2_nsteps)
      for (i in ind_lo:ind_hi) {
        tmp_ind = 2*(i - ind_lo) + 1
        curr_k_b_vec[tmp_ind] = k_b_vec[i]
        curr_par1_vec[tmp_ind] = par1_vec[i]
        curr_gdev_vec[tmp_ind] = gdev_vec[i]
      }
      for (i in 1:(ind_hi - ind_lo)) {
        tmp_ind = 2*i
        curr_k_b_vec[tmp_ind] = 0.5*(curr_k_b_vec[tmp_ind - 1] + curr_k_b_vec[tmp_ind + 1])
        curr_par1_vec[tmp_ind] = 0.5*(curr_par1_vec[tmp_ind - 1] + curr_par1_vec[tmp_ind + 1])
      }
      for (i in 1:curr_par2_nsteps) {
        if (!is.na(curr_gdev_vec[i])) { next }
        curr_k_b = data.frame(curr_k_b = curr_k_b_vec[i])
        attach(k_jam_use)
        attach(curr_k_b)
        attach(traffic_data)
        attach(inner_ccrit_use)
        attach(inner_ncyc_use)
        attach(outer_ccrit_use)
        attach(outer_ncyc_use)
        optim_obj = try(find.hyper(model = model_formula, parameters = c(curr_par1_vec[i]), k = 0.0, steps = c(par1_step), maxit = 500))
        if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
        if (optim_obj$convergence != 0) {
          par1_range = max(10.0*abs(curr_par1_vec[i]), 1000.0*par1_step)
          par1_min = curr_par1_vec[i] - par1_range
          par1_max = curr_par1_vec[i] + par1_range
          optim_obj = try(find.hyper(model = model_formula, parameters = c(curr_par1_vec[i]), k = 0.0, steps = c(par1_step), lower = c(par1_min), upper = c(par1_max),
                                     method = 'Brent', maxit = 500))
          if (class(optim_obj) == 'try-error') { optim_obj = list(convergence = 1) }
          if (optim_obj$convergence != 0) {
            cat('ERROR - The refining fits did not converge...\n')
            detach(outer_ncyc_use)
            detach(outer_ccrit_use)
            detach(inner_ncyc_use)
            detach(inner_ccrit_use)
            detach(traffic_data)
            detach(curr_k_b)
            detach(k_jam_use)
            q(save = 'no', status = 1)
          }
          if ((optim_obj$par[1] <= (par1_min + par1_step)) || (optim_obj$par[1] >= (par1_max - par1_step))) {
            cat('ERROR - The refining fits did not converge (parameter limit reached)...\n')
            detach(outer_ncyc_use)
            detach(outer_ccrit_use)
            detach(inner_ncyc_use)
            detach(inner_ccrit_use)
            detach(traffic_data)
            detach(curr_k_b)
            detach(k_jam_use)
            q(save = 'no', status = 1)
          }
        }
        detach(outer_ncyc_use)
        detach(outer_ccrit_use)
        detach(inner_ncyc_use)
        detach(inner_ccrit_use)
        detach(traffic_data)
        detach(curr_k_b)
        detach(k_jam_use)
        curr_par1_vec[i] = optim_obj$par[1]
        curr_gdev_vec[i] = optim_obj$value
      }
      curr_min_gdev = min(curr_gdev_vec)
      curr_ind_min_gdev = which(curr_gdev_vec == curr_min_gdev)
      curr_ind_min_gdev = curr_ind_min_gdev[ceiling(0.5*length(curr_ind_min_gdev))]
      par2_nsteps = curr_par2_nsteps
      k_b_vec = curr_k_b_vec
      par1_vec = curr_par1_vec
      gdev_vec = curr_gdev_vec
      ind_min_gdev = curr_ind_min_gdev
    }
    par1 = par1_vec[ind_min_gdev]
    par2 = k_b_vec[ind_min_gdev]

    # Perform the final fit
    model_obj = gamlss(V3 ~ 0 + ifelse(V2 > par2, 0.0, V2*exp(par1*V2)) + ifelse(V2 > par2, V2*log(k_jam/V2), 0.0), sigma.formula = ~ ns(V2, df = 4),
                       nu.formula = ~ ns(V2, df = 2), tau.formula = ~ 1 + V2, family = SEP3(), data = traffic_data, c.crit = outer_ccrit, n.cyc = outer_ncyc,
                       i.control = glim.control(cc = inner_ccrit, cyc = inner_ncyc))
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
    if (!all(is.finite(model_obj$nu.fv))) {
      cat('ERROR - The predicted values for "nu" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(model_obj$nu.fv <= 0.0)) {
      cat('ERROR - The predicted values for "nu" at the density values in the data include at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    if (!all(is.finite(model_obj$tau.fv))) {
      cat('ERROR - The predicted values for "tau" at the density values in the data include at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(model_obj$tau.fv <= 0.0)) {
      cat('ERROR - The predicted values for "tau" at the density values in the data include at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    traffic_data[, fitted_values_mu := model_obj$mu.fv]
    traffic_data[, fitted_values_sigma := model_obj$sigma.fv]
    traffic_data[, fitted_values_nu := model_obj$nu.fv]
    traffic_data[, fitted_values_tau := model_obj$tau.fv] },
  error = function(cond) { cat('ERROR - Failed to store the predicted values for the fitted model...\n')
                           q(save = 'no', status = 1) }
)

# Compute the normalised quantile residuals and store them in the data table. Note that the normalised quantile residuals may include some "-Inf"
# or "Inf" values for particularly bad outliers.
cat('Computing the normalised quantile residuals...\n')
tryCatch(
  { cumulative_probs_lower = pSEP3(traffic_data$V3, mu = model_obj$mu.fv, sigma = model_obj$sigma.fv, nu = model_obj$nu.fv, tau = model_obj$tau.fv)
    cumulative_probs_upper = pSEP3(traffic_data$V3, mu = model_obj$mu.fv, sigma = model_obj$sigma.fv, nu = model_obj$nu.fv, tau = model_obj$tau.fv, lower.tail = FALSE)
    traffic_data[, normalised_quantile_residuals := calculate_normalised_quantile_residuals(cumulative_probs_lower, cumulative_probs_upper)] },
  error = function(cond) { cat('ERROR - Failed to compute the normalised quantile residuals...\n')
                           q(save = 'no', status = 1) }
)

# Compute the percentiles for the fitted model at the density values in the data and store them in the data table
cat('Computing the percentiles for the fitted model at the density values in the data...\n')
tryCatch(
  { traffic_data[, percentile_m3sig := qSEP3(pNO(-3.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                             nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)]
    traffic_data[, percentile_m2sig := qSEP3(pNO(-2.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                             nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)]
    traffic_data[, percentile_m1sig := qSEP3(pNO(-1.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                             nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)]
    traffic_data[, percentile_0sig := qSEP3(0.5, mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                            nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)]
    traffic_data[, percentile_p1sig := qSEP3(pNO(1.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                             nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)]
    traffic_data[, percentile_p2sig := qSEP3(pNO(2.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                             nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)]
    traffic_data[, percentile_p3sig := qSEP3(pNO(3.0), mu = traffic_data$fitted_values_mu, sigma = traffic_data$fitted_values_sigma,
                                             nu = traffic_data$fitted_values_nu, tau = traffic_data$fitted_values_tau)] },
  error = function(cond) { cat('ERROR - Failed to compute the percentiles for the fitted model at the density values in the data...\n')
                           q(save = 'no', status = 1) }
)

# Compute a set of distributional measures for the fitted model at the density values in the data and store them in the data table
cat('Computing a set of distributional measures for the fitted model at the density values in the data...\n')
tryCatch(
  { distributional_measures = calculate_distributional_measures_for_SEP3(traffic_data$fitted_values_mu, traffic_data$fitted_values_sigma,
                                                                         traffic_data$fitted_values_nu, traffic_data$fitted_values_tau)
    traffic_data[, mean := distributional_measures$mean]
    traffic_data[, median := double(length = ntraffic_data)]
    traffic_data[, mode := distributional_measures$mode]
    traffic_data[, standard_deviation := distributional_measures$standard_deviation]
    traffic_data[, moment_skewness := distributional_measures$moment_skewness]
    traffic_data[, moment_excess_kurtosis := distributional_measures$moment_excess_kurtosis] },
  error = function(cond) { cat('ERROR - Failed to compute a set of distributional measures for the fitted model at the density values in the data...\n')
                           q(save = 'no', status = 1) }
)

# Reconstruct the fitted model over the density range from zero to "upper_density"
cat('Reconstructing the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { reconstructed_model_fit = data.table(V2 = seq(from = 0.0, to = upper_density, length.out = ngrid))
    predicted_values = predictAll(model_obj, newdata = reconstructed_model_fit, type = 'response', data = traffic_data)
    predicted_values$sigma = fix_out_of_data_curves(reconstructed_model_fit$V2, predicted_values$sigma, data_min_density, data_max_density)
    predicted_values$nu = fix_out_of_data_curves(reconstructed_model_fit$V2, predicted_values$nu, data_min_density, data_max_density)
    predicted_values$tau = fix_out_of_data_curves(reconstructed_model_fit$V2, predicted_values$tau, data_min_density, data_max_density)
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
    if (!all(is.finite(predicted_values$tau))) {
      cat('ERROR - The reconstructed fitted model for "tau" includes at least one value that is infinite...\n')
      q(save = 'no', status = 1)
    }
    if (any(predicted_values$tau <= 0.0)) {
      cat('ERROR - The reconstructed fitted model for "tau" includes at least one value that is zero or negative...\n')
      q(save = 'no', status = 1)
    }
    reconstructed_model_fit[, mu := predicted_values$mu]
    reconstructed_model_fit[, sigma := predicted_values$sigma]
    reconstructed_model_fit[, nu := predicted_values$nu]
    reconstructed_model_fit[, tau := predicted_values$tau] },
  error = function(cond) { cat('ERROR - Failed to reconstruct the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Construct percentile curves for the fitted model over the density range from zero to "upper_density"
cat('Constructing percentile curves for the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { reconstructed_model_fit[, percentile_m3sig := qSEP3(pNO(-3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                        nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)]
    reconstructed_model_fit[, percentile_m2sig := qSEP3(pNO(-2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                        nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)]
    reconstructed_model_fit[, percentile_m1sig := qSEP3(pNO(-1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                        nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)]
    reconstructed_model_fit[, percentile_0sig := qSEP3(0.5, mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                       nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)]
    reconstructed_model_fit[, percentile_p1sig := qSEP3(pNO(1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                        nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)]
    reconstructed_model_fit[, percentile_p2sig := qSEP3(pNO(2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                        nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)]
    reconstructed_model_fit[, percentile_p3sig := qSEP3(pNO(3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma,
                                                        nu = reconstructed_model_fit$nu, tau = reconstructed_model_fit$tau)] },
  error = function(cond) { cat('ERROR - Failed to construct percentile curves for the fitted model over the required density range...\n')
                           q(save = 'no', status = 1) }
)

# Construct curves of a set of distributional measures for the fitted model over the density range from zero to "upper_density"
cat('Constructing curves of a set of distributional measures for the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { distributional_measures = calculate_distributional_measures_for_SEP3(reconstructed_model_fit$mu, reconstructed_model_fit$sigma, reconstructed_model_fit$nu, reconstructed_model_fit$tau)
    reconstructed_model_fit[, mean := distributional_measures$mean]
    reconstructed_model_fit[, median := double(length = ngrid)]
    reconstructed_model_fit[, mode := distributional_measures$mode]
    reconstructed_model_fit[, standard_deviation := distributional_measures$standard_deviation]
    reconstructed_model_fit[, moment_skewness := distributional_measures$moment_skewness]
    reconstructed_model_fit[, moment_excess_kurtosis := distributional_measures$moment_excess_kurtosis] },
  error = function(cond) { cat('ERROR - Failed to construct curves of a set of distributional measures for the fitted model over the required density range...\n')
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
    npar_tau = model_obj$tau.df
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
    dvdk_0 = par1*model_obj$mu.coefficients[1]
    k_crit = NA     # Can be computed analytically - not yet implemented
    k_vmax = NA     # Can be computed analytically - not yet implemented
    q_cap = NA      # Can be computed once k_crit is available - not yet implemented
    v_max = NA      # Can be computed once k_vmax is available - not yet implemented
    v_bw = NA
    dvdk_kjam = NA
    if (model_obj$mu.coefficients[1] > 0.0) { v_ff = model_obj$mu.coefficients[1] }
    if (model_obj$mu.coefficients[2] > 0.0) {
      v_bw = model_obj$mu.coefficients[2]
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
cat('  v_ff:    ', model_obj$mu.coefficients[1], '\n')
cat('  1/k_crit:', -par1, '\n')
cat('  k_b:     ', par2, '\n')
cat('  v_bw:    ', model_obj$mu.coefficients[2], '\n')

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
        model_obj$mu.coefficients[1], '           # v_ff\n',
        -par1, '           # 1/k_crit\n',
        par2, '           # k_b\n',
        model_obj$mu.coefficients[2], '           # v_bw\n',
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
