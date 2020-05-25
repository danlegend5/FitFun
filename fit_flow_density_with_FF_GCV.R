fit_flow_density_with_FF_GCV = function(data, ngrid, upper_density, output_file1) {

# Description: This function fits a GAMLSS model to the flow-density values in "data", and it is designed to be called directly from the R script
#              "FitFun.R". The model component for the functional form of the flow-density relationship is the free-flow model (FF). The model
#              component for the noise in the flow-density relationship is defined as independent observations that follow a Gaussian distribution
#              with constant variance (GCV).
#                The input parameters "ngrid" and "upper_density" are used to define an equally spaced grid of "ngrid" density values ranging from
#              zero to "upper_density". The function employs this density grid to reconstruct the fitted model at the grid points for use in plots
#              and for estimating certain properties of the fitted model that are not directly accessible from the fitted parameter values.
#                The function creates various output files including diagnostic plots ("output_file1"; see "FitFun.R" for details). If the function
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


# Report on the GAMLSS model and the data                                                       #### FINISH ANY CLEANUPS
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
cat('  No. of flow-density measurement pairs:', ndata, '\n')
cat('  Minimum density in the data:          ', data_min_density, '\n')
cat('  Maximum density in the data:          ', data_max_density, '\n')
cat('  Minimum flow in the data:             ', data_min_flow, '\n')
cat('  Maximum flow in the data:             ', data_max_flow, '\n')
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

# Reconstruct the fitted model over the density range from zero to "upper_density"
cat('\n')
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

# Construct centile curves for the fitted model over the density range from zero to "upper_density"
cat('Constructing centile curves for the fitted model over the density range from 0 to', upper_density, '...\n')
tryCatch(
  { reconstructed_model_fit[, cent_m3sig := qNO(pNO(-3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, cent_m2sig := qNO(pNO(-2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, cent_m1sig := qNO(pNO(-1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, cent_0sig := qNO(0.5, mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, cent_p1sig := qNO(pNO(1.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, cent_p2sig := qNO(pNO(2.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)]
    reconstructed_model_fit[, cent_p3sig := qNO(pNO(3.0), mu = reconstructed_model_fit$mu, sigma = reconstructed_model_fit$sigma)] },
  error = function(cond) { cat('ERROR - Failed to construct centile curves for the fitted model over the required density range...\n')
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

# Extract information from the fit model object for the fit summary
cat('Extracting information from the fit model object for the fit summary...\n')
tryCatch(
  { npar_mu = model_obj$mu.df
    npar_sigma = model_obj$sigma.df
    npar_nu = 0
    npar_tau = 0
    npar_all = model_obj$df.fit
    gdev = model_obj$G.deviance
    aic = model_obj$aic
    bic = model_obj$sbc },
  error = function(cond) { cat('ERROR - Failed to extract information from the fit model object for the fit summary...\n')
                           q(save = 'no', status = 1) }
)

# Where possible, extract physical parameter values from the fit model object for the fit summary
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
  error = function(cond) { cat('ERROR - Failed to extract physical parameter values from the fit model object for the fit summary...\n')
                           q(save = 'no', status = 1) }
)

# Write out the fit summary file "output_file1"
cat('\n')
cat('Writing out the fit summary file:', output_file1, '\n')
write_fit_summary(output_file1, 'Flow.Density', ndata, data_min_density, data_max_density, data_min_flow, data_max_flow,
                  npar_mu, npar_sigma, npar_nu, npar_tau, npar_all, gdev, aic, bic,
                  q_0, v_ff, dvdk_0, k_crit, k_vmax, q_cap, v_max, k_jam, v_bw, dvdk_kjam,
                  curve_properties_for_mu_over_data_range, curve_properties_for_sigma_over_data_range,
                  curve_properties_for_nu_over_data_range, curve_properties_for_tau_over_data_range,
                  curve_properties_for_mu_over_full_range, curve_properties_for_sigma_over_full_range,
                  curve_properties_for_nu_over_full_range, curve_properties_for_tau_over_full_range)
cat('######################################################################################################################\n', file = output_file1, append = TRUE)
cat('# FITTED MODEL PARAMETERS (SEE THE ACCOMPANYING PAPER BY BRAMICH, MENENDEZ & AMBUHL FOR DETAILS)\n', file = output_file1, append = TRUE)
cat('# N.B: FITTED COEFFICIENTS FOR ANY NON-PARAMETRIC SMOOTHING FUNCTIONS IN THE MODEL ARE NOT REPORTED HERE\n', file = output_file1, append = TRUE)
cat('######################################################################################################################\n', file = output_file1, append = TRUE)
cat(model_obj$mu.coefficients[1], '          # v_ff\n', file = output_file1, append = TRUE)
cat(exp(model_obj$sigma.coefficients[1]), '          # sigma_con\n', file = output_file1, append = TRUE)
cat('######################################################################################################################\n', file = output_file1, append = TRUE)
cat('# FIT SUMMARY AS PROVIDED BY THE GAMLSS SOFTWARE\n', file = output_file1, append = TRUE)
cat('######################################################################################################################\n', file = output_file1, append = TRUE)
sink(file = output_file1, append = TRUE)
summary(model_obj)
sink()



#### ABOVE FULLY READ AND TESTED



q(save = 'no', status = 1)








#cat("<routes>\n", file=paste0("C:/Users/ambuehll/Networks/LIS_2/trips24h_smoothed.rou.xml"),append=F)
#write.table(trip,append = T,file=paste0("C:/Users/ambuehll/Networks/LIS_2/trips24h_smoothed.rou.xml"),col.names = F,row.names = F,quote = F)
#cat("</routes>\n", file=paste0("C:/Users/ambuehll/Networks/LIS_2/trips24h_smoothed.rou.xml"),append=T)


#The R base function write.table() can be used to export a data frame or a matrix to a file.
#A simplified format is as follow:
#write.table(x, file, append = FALSE, sep = " ", dec = ".",
#            row.names = TRUE, col.names = TRUE)
#    x: a matrix or a data frame to be written.
#    file: a character specifying the name of the result file.
#    sep: the field separator string, e.g., sep = “\t” (for tab-separated value).
#    dec: the string to be used as decimal separator. Default is “.”
#    row.names: either a logical value indicating whether the row names of x are to be written along with x, or a character vector of row names to be written.
#    col.names: either a logical value indicating whether the column names of x are to be written along with x, or a character vector of column names to
#                be written. If col.names = NA and
#                    row.names = TRUE a blank column name is added, which is the convention used for CSV files to be read by spreadsheets.


# Loading mtcars data
#data("mtcars")
# Writing mtcars data
#write.table(mtcars, file = "mtcars.txt", sep = "\t",
#            row.names = TRUE, col.names = NA)
#
#If you don’t want to write row names, use row.names = FALSE as follow:
#
#write.table(mtcars, file = "mtcars.txt", sep = "\t",
#            row.names = FALSE)




#write.table(data, output_file1)
#                           if (file.exists(output_file1)) { file.remove(output_file1) }




# Report the fit summary
cat('\n')
cat('Fit summary:\n')
cat('\n')
cat('No. of flow-density measurement pairs (Ndat):', ndata, '\n')
cat('No. of free parameters (mu):                 ', npar_mu, '\n')
cat('No. of free parameters (sigma):              ', npar_sigma, '\n')
cat('No. of free parameters (nu):                 ', npar_nu, '\n')
cat('No. of free parameters (tau):                ', npar_tau, '\n')
cat('Total no. of free parameters (Npar):         ', npar_all, '\n')
cat('Global deviance (-2 ln L):                   ', gdev, '\n')
cat('AIC (-2 ln L + 2 Npar):                      ', aic, '\n')
cat('BIC (-2 ln L + Npar ln Ndat):                ', bic, '\n')



q(save = 'no', status = 1)



#### BELOW FULLY READ AND TESTED


# Return the GAMLSS model fit object
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
return(model_obj)
}


#### ABOVE FULLY READ AND TESTED




#return(list(model_obj = model_obj, status = 1, errstr = errstr))



#library(ggplot2)
#library(gamlss.nl)
#library(gamlss.util)
#library(colorspace)
#library(ggpubr)

#plot1 = ggplot(data) + geom_point(aes(x = V2, y = V3))
#ggsave('plot1.png', plot1)

# Exclude zeros
#dataD = data[V2 > 0]
##############

#dataD[, V2_2 := V2^2]

# MODEL 1
#m0 = gamlss(V3 ~ 0 + V2 + V2_2, data = dataD, family = NO)

#summary(m0)

#model_values = data.table(V2 = seq(0.001, 1.0, 0.001))
#model_values[, V2_2 := V2^2]


#model_values[, y_pred := predict(m0, newdata = model_values, type = "response", data = dataD)]

#print(model_values)

#plot1 = ggplot() +
#        geom_point(data = dataD, aes(V2, V3), color = "red") +
#        geom_line(data = model_values, aes(V2, y_pred))
#ggsave('plot1.png', plot1)

# MODEL 1
#cat('\n')
#cat('Non-linear...\n')

#dataD[, V2_L := V2*log(V2)]

#m1 = gamlss(V3 ~ 0 + V2 + V2_L, sigma.formula = ~ pb(V2), data = dataD, family = NO)
#m1 = gamlss(V3 ~ 0 + V2 + V2_L, sigma.formula = ~ pb(V2), nu.formula = ~ pb(V2), data = dataD, family = exGAUS)

#summary(m1)

#model_values = data.table(V2 = seq(0.001, 1.0, 0.001))
#model_values[, V2_2 := V2^2]
#model_values[, V2_L := V2*log(V2)]
#model_values[, y_pred := predict(m1, newdata = model_values, type = "response", data = dataD)]

#print(model_values)
#print(dataD)

#mu_coeff = m1$mu.coefficients
#sig_coeff = m1$sigma.coefficients

#flow_dat = data.table(V3=rep(seq(1,1000,1),10),V2=rep(seq(0.07,0.87,length.out=10),each=1000))
#flow_dat[,sig := predict(m1, what = 'sigma', newdata = flow_dat, type = "response", data = dataD)]
#flow_dat[,mu := mu_coeff[1]+ mu_coeff[2]*V2*log(V2)]
#flow_dat[,flow_mod := 20.0*(1.0/sqrt(2.0*pi*sig*sig))*exp(-0.5*(((V3 - mu)/sig)^2))]


#print(flow_dat)

#plot3 = ggplot(dataD) + geom_point(aes(V2,V3),colour='grey70', size=0.5) +
#        geom_line(data = model_values, aes(V2, y_pred), size = 1.0) +
#        geom_path(data = flow_dat, aes(V2 - flow_mod,V3,group=V2), colour='red') +
#        theme_pubr(base_size = 16, border = TRUE) + scale_y_continuous(expand = c(0,0)) +
#        xlab('Occupancy') + ylab('Flow (veh/hour)') + ggtitle('uk : london : primary : CNTR_N01_523b1') +
#        theme(plot.title = element_text(hjust = 0.5)) +
#        geom_vline(xintercept = seq(0.07,0.87,length.out=10), linetype = 'dashed', color = 'red')

#ggsave('plot3.png', plot3, width = 6.0, height = 4.0, scale = 2.0)


#q()


#plotSimpleGamlss(V3, exp(V2_L), model = m1, data = dataD, x.val = seq(0.0, 1.0, 0.2), val = 5, N = 1000, ylim = c(0.0, 1000.0))

#q()

#plot2 = ggplot() +
#        geom_point(data = dataD, aes(V2, V3), color = "red") +
#        geom_line(data = model_values, aes(V2, y_pred))
#ggsave('plot2.png', plot2)

#q()

# MODEL 2

#cat('\n')
#cat('Non-linear...\n')

#m2 = nlgamlss(y = V3, mu.formula = ~ V2 + V2_2, mu.start = c(0.0, 1000.0, -1000.0), sigma.start = 1.0, family = NO, data = data)
#summary(m2)


#q()


# fit

#res <- gamlss(y~cs(data$V2, df = 5), sigma.fo = ~cs(data$V2, df = 5), family = NO, data = data$V3)

#summary(res)

# predict on grid

#pred <- data.table(occ=seq(0,0.8,0.001))

#pred <- pred[,y_pred:=predict(res,newdata = pred,type="response")]



# fit on x (orignal data)

#tt[,y_hat:=fitted(res)]

#tt[,residuals:=residuals(res)]




#ggplot(tt)+geom_point(aes(occ,flow),color="red")+

  # geom_line(aes(occ,y_hat))+

#  geom_line(data=pred,aes(occ,y_pred),colour="blue")
