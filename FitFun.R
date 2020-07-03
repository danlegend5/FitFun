# FitFun: Fitting Fundamental Diagrams Of Road Traffic
#
# Description: ???This script???  Cite the paper and book "Flexible Regression and Smoothing: Using GAMLSS in R" by Stasinopoulos et al.       #### FINISH
#                The specific procedures for fitting each individual model that is available in this script are separated out into a single R
#              module per model. These modules have names of the form "fit_flow_density_with_<functional_form_model>_<noise_model>.R" or
#              "fit_speed_density_with_<functional_form_model>_<noise_model>.R", and they each contain a single function of the same name (see
#              below under "Command-Line Arguments" for the definitions of "functional_form_model" and "noise_model"). Any parameters that are
#              specific to a procedure for fitting a particular model are defined within the relevant module itself, and their default values
#              can be modified if necessary.
#
# Important Notes On Model Selection Criteria:
#
#     Very often in scientific modelling, one faces the task of selecting an optimal (or best) model for a data sample from a set of candidate
#   models (multiple working hypotheses). In this context, "optimality" refers both to the Principle of Parsimony, in that the best model should
#   constitute the simplest model that provides a good fit to the data without under- or over-fitting, and to appropriate/relevant model
#   performance measure(s). Model estimation via maximum likelihood assumes a uniform prior probability density function on the model parameters.
#   Consequently, as parameters are added to a model, the maximum likelihood always increases, rendering it useless for the purpose of model
#   selection between models with different dimensionality. Information criteria are used as an alternative for evaluating models with different
#   numbers of parameters. Various information criteria have been developed from distinct statistical view-points as implementations of the
#   Principle of Parsimony and each one may be used to automatically select a parsimonious model from a set of candidate models. They may be
#   applied regardless of whether the models under consideration are nested or non-nested.
#     The GAMLSS software, on which this script is based, computes both the Akaike information criterion (AIC; Akaike 1974, IEEE Transactions on
#   Automatic Control, 19, 716) and the Bayesian information criterion (BIC; Schwarz 1978, The Annals of Statistics, 6, 461) for each GAMLSS
#   model fit. Both of these information criteria apply to linear and non-linear models estimated via maximum likelihood. The formulae for these
#   information criteria are given by:
#
#   AIC = -2*ln(Lmax) + 2*Npar
#
#   BIC = -2*ln(Lmax) + Npar*ln(Ndat)
#
#   where Lmax is the value of the likelihood function L at the vector of maximum likelihood estimators (MLEs) for the model parameters, Npar is
#   the number of free parameters in the model, and Ndat is the number of data values. Model selection with the AIC or BIC is performed by
#   minimising -2*ln(L) for each model, and then minimising AIC or BIC, respectively, over the full set of models under consideration.
#     The AIC is derived as an asymptotic approximation to the Kullback-Leibler divergence (Kullback & Leibler 1951, The Annals of Mathematical
#   Statistics, 22, 79) which measures the distance of a candidate model from the true underlying model under the (reasonable) assumption that
#   the true model is of infinite dimension. In this case the AIC provides an asymptotically efficient selection of a finite dimensional
#   approximating model. However, if the true model is finite dimensional, then the AIC does not provide a consistent model selection. A model
#   selection criterion is said to be consistent if it selects with high probability the true model from the set of candidate models whenever
#   the true model is represented in the set of candidate models. The aim of the AIC is to evaluate models based on their prediction accuracy.
#     It is important to be aware that the AIC suffers from a large negative bias for small samples, or when the number of model parameters is a
#   non-negligible fraction of the number of data values (e.g. for Ndat/Npar < 40; Hurvich & Tsai 1989, Biometrika, 76, 297). The bias in the AIC
#   can be corrected at the cost of its general applicability, since the formula for the bias correction depends on the statistical model being
#   fitted. For example, Sugiura (1978, Communications in Statistics - Theory and Methods, 7, 13) derived a bias-corrected version of the AIC for
#   Gaussian linear regression problems that is asymptotically the same as the AIC for Ndat >> Npar.
#     Takeuchi (1976, Mathematical Sciences, 153, 12) generalised the AIC with a more complicated formula to create the Takeuchi information
#   criterion (TIC). Subsequently, Konishi & Kitagawa (1996, Biometrika, 83, 875) derived a further generalisation of the AIC and TIC, called the
#   generalised information criterion (GIC), that can also be applied to model selection for models with parameters estimated by maximum penalised
#   likelihood. The GAMLSS software only implements the AIC from the available AIC-like information criteria. For GAMLSS models that employ
#   non-parametric smoothing terms estimated by maximum penalised likelihood, the AIC is computed using the trace of the smoother matrix as an
#   estimate of the (effective) number of free parameters used in the fit.
#     An alternative approach to model selection is a Bayesian approach where the model with the largest Bayesian posterior probability is chosen.
#   The BIC is derived by approximating the posterior probability of each model. The BIC generally includes a heavier penalty than the AIC (and
#   the small-sample bias-corrected AIC) for more complicated models (e.g. in the regime Npar < 20 for Ndat > 50), therefore favouring models with
#   fewer parameters than those favoured by the AIC. Konishi, Ando & Imoto (2004, Biometrika, 91, 27) performed a deeper Bayesian analysis to
#   derive an improved BIC, along with a version that applies to model selection for models with parameters estimated by maximum penalised
#   likelihood. The BIC and the improved BIC are consistent model selection criteria. Again, the GAMLSS software only implements the BIC from the
#   available BIC-like information criteria, and for GAMLSS models that employ non-parametric smoothing terms estimated by maximum penalised
#   likelihood, the trace of the smoother matrix is adopted as an estimate of the effective number of free parameters.
#
# Usage:
#
#   To run this script, issue the following command:
#
#   <path_to_binary>/Rscript --vanilla <path_to_script>/FitFun.R <path_to_modules> <data_file> <output_dir> <overwrite> <fd_type>
#                                      <functional_form_model> <noise_model> <ngrid> <upper_density> <plot_format>
#
#   where "path_to_binary" is the full directory path to the "Rscript" binary (usually "/usr/bin") and "path_to_script" is the full directory
#   path to where the script "FitFun.R" is stored. The definition of each of the command-line arguments "path_to_modules", "data_file",
#   "output_dir", "overwrite", "fd_type", "functional_form_model", "noise_model", "ngrid", "upper_density", and "plot_format" can be found below.
#   If "path_to_binary" is in the user's path, then it can be dropped. Also, if the script is being run from the directory where it resides, then
#   "path_to_script" can also be dropped. In this case, the command reduces to:
#
#   Rscript --vanilla FitFun.R <path_to_modules> <data_file> <output_dir> <overwrite> <fd_type> <functional_form_model> <noise_model> <ngrid>
#                              <upper_density> <plot_format>
#
# Command-Line Arguments:
#
#   path_to_modules - STRING - The full directory path to where the R modules required by this script are stored.
#   data_file - STRING - File name of the input data file. This argument can be supplied with or without a full directory path. If this argument
#                        is supplied without a full directory path, then the script will look for the input data file in the directory from where
#                        the script is run. See below for the required format of the data file.
#   output_dir - STRING - The output directory where the output files are to be written. This argument can be supplied with or without a full
#                         directory path. If this argument is supplied without a full directory path, then the script will assume that the output
#                         directory is inside the directory from where the script is run. The following output files will be written to the
#                         output directory:
#
#                         Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt
#                         Fit.Curves.<fd_type>.<functional_form_model>.<noise_model>.txt
#                         Fit.Predictions.<fd_type>.<functional_form_model>.<noise_model>.txt
#                         Plot.Of.Fitted.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Residuals.From.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Percentiles.And.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Normalised.Quantile.Residuals.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Fitted.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Residuals.From.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Percentiles.And.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Normalised.Quantile.Residuals.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Normalised.Quantile.Residuals.Versus.Mu.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Normalised.Quantile.Residuals.Versus.Time.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Detrended.Normal.QQ.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                         Plot.Of.Slotted.ACF.For.Normalised.Quantile.Residuals.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#
#                         See below for a description of each of the output files.
#   overwrite - STRING - If this argument is set to 'yes', then the script will overwrite any of the output files that already exist in the output
#                        directory "output_dir". If this argument is set to 'no', then the script will stop without doing anything if it finds
#                        that any of the output files already exist in the output directory "output_dir". If this argument is set to any other
#                        string, then the script will fail.
#   fd_type - STRING - The form of the empirical FD data as a string in the format "<dependent_variable>.<independent_variable>". The acceptable
#                      values for this argument are 'Flow.Density' and 'Speed.Density'. If this argument is set to any other string, then the
#                      script will fail.
#   functional_form_model - STRING - ???model component for the functional form???                                                              #### FINISH
#   noise_model - STRING - ???model component for the noise???                                                                                  #### FINISH
#   ngrid - INTEGER - The script will define an equally spaced grid of "ngrid" density values for the purpose of reconstructing the fitted model.
#                     This argument must be a positive number greater than or equal to "11". It is recommended to set this argument to at least
#                     "10001".
#   upper_density - FLOAT - The script will define the grid of equally spaced density values to cover the range from zero to "upper_density" for
#                           the purpose of reconstructing the fitted model. This argument must be a positive number, and it must also be greater
#                           than or equal to the maximum value of the independent variable (i.e. density or occupancy) in the data.
#   plot_format - STRING - The file format for the plot files that the script produces. The acceptable values for this argument are 'ps', 'pdf',
#                          'png', and 'none'. If this argument is set to any other string, then the script will fail. In the case that this
#                          argument is set to 'none', then the script will not produce any plot files.
#
# Input Data File:
#
#   The input data file "data_file" should be a non-empty ASCII text file with exactly three columns separated by white space of any length. The
#   data file should not contain any header lines, or any lines that are not data lines, and it should have at least 5 data lines. The required
#   columns are:
#
#   Column 1 - INTEGER/FLOAT - A time stamp corresponding to the start/midpoint/end of the time interval over which the traffic measurements were
#                              taken. The units of time, and the point within the time interval to which the time stamp corresponds, are not
#                              important, so long as they are consistent over all of the traffic measurement intervals, since this allows reliable
#                              time differences to be computed. Currently, the data in this column are only used for the computation of the slotted
#                              auto-correlation function for the normalised quantile residuals.
#   Column 2 - INTEGER/FLOAT - The measured value of the independent variable (i.e. density or occupancy) in the corresponding time interval. All
#                              values in this column must be positive.
#   Column 3 - INTEGER/FLOAT - The measured value of the dependent variable (i.e. flow or speed) in the corresponding time interval. All values in
#                              this column must be non-negative, and there must be at least one value that is non-zero.
#
# Output Files:
#
#   Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt - This output text file provides a set of summary information for the fitted
#                                                                     model. The contents of the file are fully documented within the file itself.
#   Fit.Curves.<fd_type>.<functional_form_model>.<noise_model>.txt - This output text file provides the reconstructed fitted model, along with
#                                                                    percentile curves, for the equally spaced grid of "ngrid" density values
#                                                                    covering the range from zero to "upper_density". A header line provides the
#                                                                    column descriptions.
#   Fit.Predictions.<fd_type>.<functional_form_model>.<noise_model>.txt - This output text file provides the predicted values for the model,
#                                                                         along with the normalised quantile residuals, at the density values in
#                                                                         the data. A header line provides the column descriptions. For more
#                                                                         information on what a normalised quantile residual is, please see
#                                                                         Chapter 12 in the book "Flexible Regression and Smoothing: Using GAMLSS
#                                                                         in R" by Stasinopoulos et al.
#   Plot.Of.Fitted.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the flow or speed data versus density (red points). The fitted model
#                                                                    component for "mu" over the density range from zero to the maximum observed
#                                                                    density is plotted as a black curve.
#   Plot.Of.Residuals.From.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the residuals of the flow or speed data from the fitted model
#                                                                    component for "mu" versus density (red points). The plot density range is
#                                                                    from zero to the maximum observed density.
#   Plot.Of.Percentiles.And.Mu.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the flow or speed data versus density (red points). The fitted model
#                                                                    component for "mu" over the density range from zero to the maximum observed
#                                                                    density is plotted as a black curve. The percentile ranges corresponding to
#                                                                    0.135 to 2.28, 2.28 to 15.87, and 15.87 to 50.0 are plotted as lighter to
#                                                                    darker grey regions, respectively, below the median curve (dashed blue curve).
#                                                                    The median curve may not be visible if it is coincident with the "mu" curve.
#                                                                    The percentile ranges corresponding to 50.0 to 84.13, 84.13 to 97.72, and 97.72
#                                                                    to 99.865 are plotted as darker to lighter grey regions, respectively, above
#                                                                    the median curve. The percentile boundaries correspond to -3*sigma, -2*sigma,
#                                                                    -1*sigma, the median, 1*sigma, 2*sigma, and 3*sigma in a Normal distibution.
#   Plot.Of.Normalised.Quantile.Residuals.For.Data.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the normalised quantile residuals for the flow or speed data versus
#                                                                    density (red points). The plot density range is from zero to the maximum
#                                                                    observed density. The grey regions are percentile ranges as described for
#                                                                    previous plots while the median line is coincident with the horizontal dotted
#                                                                    line.
#   Plot.Of.Fitted.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the flow or speed data versus density (red points). The fitted model
#                                                                    component for "mu" over the density range from zero to "upper_density" is
#                                                                    plotted as a black curve.
#   Plot.Of.Residuals.From.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the residuals of the flow or speed data from the fitted model
#                                                                    component for "mu" versus density (red points). The plot density range is
#                                                                    from zero to "upper_density".
#   Plot.Of.Percentiles.And.Mu.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the flow or speed data versus density (red points). The fitted model
#                                                                    component for "mu" over the density range from zero to "upper_density" is
#                                                                    plotted as a black curve. The percentile ranges corresponding to 0.135 to 2.28,
#                                                                    2.28 to 15.87, and 15.87 to 50.0 are plotted as lighter to darker grey regions,
#                                                                    respectively, below the median curve (dashed blue curve). The median curve may
#                                                                    not be visible if it is coincident with the "mu" curve. The percentile ranges
#                                                                    corresponding to 50.0 to 84.13, 84.13 to 97.72, and 97.72 to 99.865 are plotted
#                                                                    as darker to lighter grey regions, respectively, above the median curve. The
#                                                                    percentile boundaries correspond to -3*sigma, -2*sigma, -1*sigma, the median,
#                                                                    1*sigma, 2*sigma, and 3*sigma in a Normal distibution.
#   Plot.Of.Normalised.Quantile.Residuals.For.Full.Density.Range.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the normalised quantile residuals for the flow or speed data versus
#                                                                    density (red points). The plot density range is from zero to "upper_density".
#                                                                    The grey regions are percentile ranges as described for previous plots while
#                                                                    the median line is coincident with the horizontal dotted line.
#   Plot.Of.Normalised.Quantile.Residuals.Versus.Mu.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the normalised quantile residuals for the flow or speed data versus
#                                                                    the fitted values for "mu" (red points). The grey regions are percentile
#                                                                    ranges as described for previous plots while the median line is coincident
#                                                                    with the horizontal dotted line.
#   Plot.Of.Normalised.Quantile.Residuals.Versus.Time.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the normalised quantile residuals for the flow or speed data versus
#                                                                    time (red points). The grey regions are percentile ranges as described for
#                                                                    previous plots while the median line is coincident with the horizontal dotted
#                                                                    line.
#   Plot.Of.Detrended.Normal.QQ.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Given a relevant set of theoretical quantiles from a standard Normal distribution,
#                                                                    this output file displays a plot of the deviation of the normalised quantile
#                                                                    residuals from the theoretical quantiles (units of sigma) versus the theoretical
#                                                                    quantiles (units of sigma; red points). The 95% confidence interval is plotted
#                                                                    as a pair of dashed curves. This sort of plot is referred to as a detrended
#                                                                    Normal quantile-quantile plot (or worm plot). For more details, see van Buuren
#                                                                    & Fredriks (2001, Statistics in Medicine, 20, 1259).
#   Plot.Of.Slotted.ACF.For.Normalised.Quantile.Residuals.<fd_type>.<functional_form_model>.<noise_model>.<plot_format>
#                                                                  - Plot of the slotted auto-correlation function (slotted ACF) versus time lag
#                                                                    (light grey bars) for the normalised quantile residuals. Uncertainties on the
#                                                                    slotted ACF values are plotted as error bars. For more details, see Edelson &
#                                                                    Krolik (1988, Astrophysical Journal, 333, 646).
#
# Requirements:
#
#   R (Version >= 3.6.1)
#   R Package: data.table (Version >= 1.12.8)
#   R Package: gamlss (Version >= 5.1-6)
#   R Package: ggplot2 (Version >= 3.2.1) [Not required if "plot_format" is set to 'none']
#   R Package: ggpubr (Version >= 0.2.5)  [Not required if "plot_format" is set to 'none']
#
# Authors:
#
#   Dan Bramich (dan.bramich@hotmail.co.uk)
#   Lukas Ambuhl (lukas.ambuehl@ivt.baug.ethz.ch)


# Ingest the command-line arguments
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('FitFun: Fitting Fundamental Diagrams Of Road Traffic\n')
cat('\n')
cat('Dan Bramich & Lukas Ambuhl\n')
cat('\n')
cat('>-----------------------------------------------------------------------------<\n')
cat('\n')
cat('Ingesting the command-line arguments...\n')
args = commandArgs(trailingOnly = TRUE)
nargs = length(args)
if (nargs < 10) {
  cat('ERROR - Too few command-line arguments specified...\n')
  q(save = 'no', status = 1)
} else if (nargs > 10) {
  cat('ERROR - Too many command-line arguments specified...\n')
  q(save = 'no', status = 1)
} else {
  path_to_modules = args[1]
  data_file = args[2]
  output_dir = args[3]
  overwrite = args[4]
  fd_type = args[5]
  functional_form_model = args[6]
  noise_model = args[7]
  ngrid = as.integer(args[8])
  upper_density = as.double(args[9])
  plot_format = args[10]
}

# Perform checks on the command-line arguments and report their values
cat('\n')
cat('Path to R modules:                                  ', path_to_modules, '\n')
if (!dir.exists(path_to_modules)) {
  cat('ERROR - The directory path to the required R modules does not exist...\n')
  q(save = 'no', status = 1)
}
cat('Input data file:                                    ', data_file, '\n')
if (!file.exists(data_file)) {
  cat('ERROR - The input data file does not exist...\n')
  q(save = 'no', status = 1)
}
cat('Output directory:                                   ', output_dir, '\n')
if ((overwrite != 'yes') && (overwrite != 'no')) {
  cat('ERROR - The command-line argument "overwrite" does not have an acceptable value...\n')
  q(save = 'no', status = 1)
}
cat('Overwrite previous files:                           ', overwrite, '\n')
acceptable_values = c('Flow.Density', 'Speed.Density')
if (!is.element(fd_type, acceptable_values)) {
  cat('ERROR - The command-line argument "fd_type" does not have an acceptable value...\n')
  q(save = 'no', status = 1)
}
cat('Fundamental diagram type:                           ', fd_type, '\n')
acceptable_values = c('FF',         'GS1935',     'GS1935kjf', 'GB1959',     'GB1959kjf',  'ED1961',     'ED1961kjf',  'UW1961A',   'UW1961B',
                      'UW1961Bkjf', 'NW1961',     'NW1961kjf', 'GZ1961A',    'GZ1961Akjf', 'GZ1961B',    'GZ1961Bkjf', 'GZ1961C',   'GZ1961Ckjf',
                      'GZ1961D',    'GZ1961Dkjf', 'GZ1961E',   'GZ1961Ekjf', 'GZ1961F',    'GZ1961G',    'GZ1961Gkjf', 'GZ1961H',   'GZ1961Hkjf',
                      'DK1966A',    'DK1966Akjf', 'DK1966B',   'DK1966Bkjf', 'DK1966C',    'DK1966Ckjf', 'MJ1971',     'MJ1971kjf', 'BM1977',
                      'VA1995',     'VA1995kjf',  'BD1995',    'DC1995A',    'DC1995Akjf', 'DC2012B',    'DC2012Bkjf', 'SN2014')                          #### FINISH
if (!is.element(functional_form_model, acceptable_values)) {
  cat('ERROR - The command-line argument "functional_form_model" does not have an acceptable value...\n')
  q(save = 'no', status = 1)
}
cat('Functional form model:                              ', functional_form_model, '\n')
acceptable_values = c('GCV', 'gaussian_quadratic_var', 'gaussian_bspline_var')                                                             #### FINISH
if (!is.element(noise_model, acceptable_values)) {
  cat('ERROR - The command-line argument "noise_model" does not have an acceptable value...\n')
  q(save = 'no', status = 1)
}
cat('Noise model:                                        ', noise_model, '\n')
if (is.na(ngrid)) {
  cat('ERROR - The command-line argument "ngrid" is not an integer...\n')
  q(save = 'no', status = 1)
}
if (ngrid < 11) {
  cat('ERROR - The command-line argument "ngrid" is an integer with a value that is less than "11"...\n')
  q(save = 'no', status = 1)
}
cat('No. of density grid points for model reconstruction:', ngrid, '\n')
if (is.na(upper_density)) {
  cat('ERROR - The command-line argument "upper_density" is not a number...\n')
  q(save = 'no', status = 1)
}
if (upper_density <= 0.0) {
  cat('ERROR - The command-line argument "upper_density" is a number that is zero or negative...\n')
  q(save = 'no', status = 1)
}
cat('Upper density for model reconstruction:             ', upper_density, '\n')
acceptable_values = c('ps', 'pdf', 'png', 'none')
if (!is.element(plot_format, acceptable_values)) {
  cat('ERROR - The command-line argument "plot_format" does not have an acceptable value...\n')
  q(save = 'no', status = 1)
}
cat('Plot file format:                                   ', plot_format, '\n')

# Load the required R libraries and functions
cat('\n')
cat('Loading the "data.table" R library...\n')
tryCatch(
  { library(data.table) },
  error = function(cond) { cat('ERROR - Failed to load the "data.table" R library...\n')
                           q(save = 'no', status = 1) }
)
cat('Loading the "gamlss" R library...\n')
cat('----------------------------------------------------------------\n')
tryCatch(
  { library(gamlss) },
  error = function(cond) { cat('ERROR - Failed to load the "gamlss" R library...\n') 
                           cat('----------------------------------------------------------------\n')
                           q(save = 'no', status = 1) }
)
cat('----------------------------------------------------------------\n')
if (plot_format != 'none') {
  cat('Loading the "ggplot2" R library...\n')
  tryCatch(
    { library(ggplot2) },
    error = function(cond) { cat('ERROR - Failed to load the "ggplot2" R library...\n')
                             q(save = 'no', status = 1) }
  )
  cat('Loading the "ggpubr" R library...\n')
  cat('----------------------------------------------------------------\n')
  tryCatch(
    { library(ggpubr) },
    error = function(cond) { cat('ERROR - Failed to load the "ggpubr" R library...\n')
                             cat('----------------------------------------------------------------\n')
                             q(save = 'no', status = 1) }
  )
  cat('----------------------------------------------------------------\n')
}
cat('Loading some R functions specific to this script...\n')
tryCatch(
  { source(file.path(path_to_modules, 'fitfun_functions.R')) },
  error = function(cond) { cat('ERROR - Failed to load the R functions...\n')
                           q(save = 'no', status = 1) }
)

# Define the names of the output files
tmpstr1 = paste0(fd_type, '.', functional_form_model, '.', noise_model)
output_files = character(length = 15)
output_files[1] = file.path(output_dir, paste0('Fit.Summary.', tmpstr1, '.txt'))
output_files[2] = file.path(output_dir, paste0('Fit.Curves.', tmpstr1, '.txt'))
output_files[3] = file.path(output_dir, paste0('Fit.Predictions.', tmpstr1, '.txt'))
if (plot_format == 'none') {
  output_files = output_files[1:3]
} else {
  tmpstr2 = paste0(tmpstr1, '.', plot_format)
  output_files[4] = file.path(output_dir, paste0('Plot.Of.Fitted.Mu.For.Data.Density.Range.', tmpstr2))
  output_files[5] = file.path(output_dir, paste0('Plot.Of.Residuals.From.Mu.For.Data.Density.Range.', tmpstr2))
  output_files[6] = file.path(output_dir, paste0('Plot.Of.Percentiles.And.Mu.For.Data.Density.Range.', tmpstr2))
  output_files[7] = file.path(output_dir, paste0('Plot.Of.Normalised.Quantile.Residuals.For.Data.Density.Range.', tmpstr2))
  output_files[8] = file.path(output_dir, paste0('Plot.Of.Fitted.Mu.For.Full.Density.Range.', tmpstr2))
  output_files[9] = file.path(output_dir, paste0('Plot.Of.Residuals.From.Mu.For.Full.Density.Range.', tmpstr2))
  output_files[10] = file.path(output_dir, paste0('Plot.Of.Percentiles.And.Mu.For.Full.Density.Range.', tmpstr2))
  output_files[11] = file.path(output_dir, paste0('Plot.Of.Normalised.Quantile.Residuals.For.Full.Density.Range.', tmpstr2))
  output_files[12] = file.path(output_dir, paste0('Plot.Of.Normalised.Quantile.Residuals.Versus.Mu.', tmpstr2))
  output_files[13] = file.path(output_dir, paste0('Plot.Of.Normalised.Quantile.Residuals.Versus.Time.', tmpstr2))
  output_files[14] = file.path(output_dir, paste0('Plot.Of.Detrended.Normal.QQ.', tmpstr2))
  output_files[15] = file.path(output_dir, paste0('Plot.Of.Slotted.ACF.For.Normalised.Quantile.Residuals.', tmpstr2))
}

# If the output directory "output_dir" already exists
if (dir.exists(output_dir)) {

  # Report
  cat('\n')
  cat('The output directory already exists:', output_dir, '\n')
  cat('Checking if any of the output files already exist...\n')

  # For each output file
  for (outfile in output_files) {

    # If the current output file already exists
    if (file.exists(outfile)) {

      # If the command-line argument "overwrite" is set to 'yes'
      if (overwrite == 'yes') {

        # Remove the current output file
        cat('Removing the output file:', outfile, '\n')
        tryCatch(
          { file.remove(outfile) },
          error = function(cond) { cat('ERROR - Failed to remove the output file...\n')
                                   q(save = 'no', status = 1) }
        )

      # If the command-line argument "overwrite" is set to 'no'
      } else {

        # Stop the script without doing anything
        cat('ERROR - The following output file already exists:', outfile, '\n')
        q(save = 'no', status = 1)
      }
    }
  }

# If the output directory "output_dir" does not already exist
} else {

  # Create the output directory "output_dir"
  cat('\n')
  cat('Creating the output directory:', output_dir, '\n')
  tryCatch(
    { dir.create(output_dir) },
    error = function(cond) { cat('ERROR - Failed to create the output directory...\n')
                             q(save = 'no', status = 1) }
  )
}

# Read in the data file "data_file"
cat('\n')
cat('Reading in the data file:', data_file, '\n')
tryCatch(
  { traffic_data = fread(data_file, header = FALSE) },
  error = function(cond) { cat('ERROR - Failed to read in the data file...\n')
                           q(save = 'no', status = 1) }
)

# Perform basic checks on the data that have been read in
ntraffic_data = nrow(traffic_data)
if (ntraffic_data == 0) {
  cat('ERROR - The data file is empty...\n')
  q(save = 'no', status = 1)
}
if (ncol(traffic_data) != 3) {
  cat('ERROR - The data file does not have exactly three columns...\n')
  q(save = 'no', status = 1)
}
if (!is.numeric(traffic_data$V1)) {
  cat('ERROR - The data in the first column of the data file are not numeric...\n')
  q(save = 'no', status = 1)
}
if (!is.numeric(traffic_data$V2)) {
  cat('ERROR - The data in the second column of the data file are not numeric...\n')
  q(save = 'no', status = 1)
}
if (!is.numeric(traffic_data$V3)) {
  cat('ERROR - The data in the third column of the data file are not numeric...\n')
  q(save = 'no', status = 1)
}
if (ntraffic_data < 5) {
  cat('ERROR - The data file has less than 5 data lines...\n')
  q(save = 'no', status = 1)
}
if (!all(is.finite(traffic_data$V1))) {
  cat('ERROR - At least one data value in the first column is infinite...\n')
  q(save = 'no', status = 1)
}
if (!all(is.finite(traffic_data$V2))) {
  cat('ERROR - At least one data value in the second column is infinite...\n')
  q(save = 'no', status = 1)
}
if (!all(is.finite(traffic_data$V3))) {
  cat('ERROR - At least one data value in the third column is infinite...\n')
  q(save = 'no', status = 1)
}
if (any(traffic_data$V2 <= 0.0)) {
  cat('ERROR - At least one data value in the second column is zero or negative...\n')
  q(save = 'no', status = 1)
}
if (any(traffic_data$V3 < 0.0)) {
  cat('ERROR - At least one data value in the third column is negative...\n')
  q(save = 'no', status = 1)
}
if (all(traffic_data$V3 == 0.0)) {
  cat('ERROR - All of the data values in the third column are zero...\n')
  q(save = 'no', status = 1)
}
cat('No. of lines read in:    ', ntraffic_data, '\n')

# Sort the data into increasing time order
cat('Sorting the data by time...\n')
tryCatch(
  { setorder(traffic_data, V1) },
  error = function(cond) { cat('ERROR - Failed to sort the data...\n')
                           q(save = 'no', status = 1) }
)

# Make a further check on "upper_density"
if (upper_density < max(traffic_data$V2)) {
  cat('ERROR - The command-line argument "upper_density" is less than the maximum value of the density data...\n')
  q(save = 'no', status = 1)
}


################################################################################################################################################
# Flow-density fundamental diagrams

# If the form of the FD relationship is flow-density
cat('\n')
cat('Calling the required R module for performing the fit...\n')
if (fd_type == 'Flow.Density') {

  # Load the required R module for performing the fit
  tryCatch(
    { source(file.path(path_to_modules, paste0('fit_flow_density_with_', functional_form_model, '_', noise_model, '.R'))) },
    error = function(cond) { cat('ERROR - Failed to load the required R module...\n')
                             q(save = 'no', status = 1) }
  )

  # Fit the chosen GAMLSS model to the data
  tryCatch(
    { fit_function = get(paste0('fit_flow_density_with_', functional_form_model, '_', noise_model))
      model_obj = fit_function(traffic_data, ngrid, upper_density, output_files) },
    error = function(cond) { cat('ERROR - Failed to fit the GAMLSS model for unknown reasons...\n')
                             remove_file_list(output_files)
                             q(save = 'no', status = 1) }
  )


################################################################################################################################################
# Speed-density fundamental diagrams

# If the form of the FD relationship is speed-density
} else if (fd_type == 'Speed.Density') {

  # Load the required R module for performing the fit
  tryCatch(
    { source(file.path(path_to_modules, paste0('fit_speed_density_with_', functional_form_model, '_', noise_model, '.R'))) },
    error = function(cond) { cat('ERROR - Failed to load the required R module...\n')
                             q(save = 'no', status = 1) }
  )

  # Fit the chosen GAMLSS model to the data
  tryCatch(
    { fit_function = get(paste0('fit_speed_density_with_', functional_form_model, '_', noise_model))
      model_obj = fit_function(traffic_data, ngrid, upper_density, output_files) },
    error = function(cond) { cat('ERROR - Failed to fit the GAMLSS model for unknown reasons...\n')
                             remove_file_list(output_files)
                             q(save = 'no', status = 1) }
  )
}

# Finish
cat('\n')
cat('Finished!\n')
