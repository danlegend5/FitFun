# FitFun: Fitting Fundamental Diagrams Of Road Traffic
#
# Description: ???This script???                          #### FINISH
#                The specific procedures for fitting each individual model that is available in this script are separated out into a single R
#              module per model. These modules have names of the form "fit_flow_density_with_<functional_form_model>_<noise_model>.R" or
#              "fit_speed_density_with_<functional_form_model>_<noise_model>.R", and they each contain a single function of the same name (see
#              below under "Command-Line Arguments" for the definitions of <functional_form_model> and <noise_model>). Any parameters that are
#              specific to a procedure for fitting a particular model are defined within the relevant module itself, and their default values
#              can be modified if necessary.
#
# Usage:
#
#   To run this script, issue the following command:
#
#   <path_to_binary>/Rscript --vanilla <path_to_script>/FitFun.R <path_to_modules> <data_file> <output_dir> <overwrite> <fd_type>
#                                      <functional_form_model> <noise_model> <ngrid> <upper_density>
#
#   where <path_to_binary> is the full directory path to the "Rscript" binary (usually "/usr/bin") and <path_to_script> is the full directory
#   path to where the script "FitFun.R" is stored. The definition of each of the command-line arguments <path_to_modules>, <data_file>,
#   <output_dir>, <overwrite>, <fd_type>, <functional_form_model>, <noise_model>, <ngrid>, and <upper_density> can be found below. If
#   <path_to_binary> is in the user's path, then it can be dropped. Also, if the script is being run from the directory where it resides, then
#   <path_to_script> can also be dropped. In this case, the command reduces to:
#
#   Rscript --vanilla FitFun.R <path_to_modules> <data_file> <output_dir> <overwrite> <fd_type> <functional_form_model> <noise_model> <ngrid>
#                              <upper_density>
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
#                         ???               #### FINISH
#                         ???               #### FINISH
#
#                         See below for a description of each of the output files.
#   overwrite - STRING - If this argument is set to 'yes', then the script will overwrite any of the output files that already exist in the output
#                        directory "output_dir". If this argument is set to 'no', then the script will stop without doing anything if it finds
#                        that any of the output files already exist in the output directory "output_dir". If this argument is set to any other
#                        string, then the script will fail.
#   fd_type - STRING - The form of the empirical FD data as a string in the format "<dependent_variable>.<independent_variable>". The acceptable
#                      values for this argument are 'Flow.Density' and 'Speed.Density'. If this argument is set to any other string, then the
#                      script will fail.
#   functional_form_model - STRING - ???model component for the functional form???     #### FINISH
#   noise_model - STRING - ???model component for the noise???                         #### FINISH
#   ngrid - INTEGER - The script will define an equally spaced grid of "ngrid" density values for the purpose of reconstructing the fitted model.
#                     This argument must be a positive number greater than or equal to "11". It is recommended to set this argument to at least
#                     "10001".
#   upper_density - FLOAT - The script will define the grid of equally spaced density values to cover the range from zero to "upper_density" for
#                           the purpose of reconstructing the fitted model. This argument must be a positive number, and it must also be greater
#                           than or equal to the maximum value of the independent variable (i.e. density or occupancy) in the data.
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
#                              time differences to be computed. (??? WOULD ONLY BE USED FOR NOISE MODEL WITH COVARIANCE ???)                            #### FINISH
#   Column 2 - INTEGER/FLOAT - The measured value of the independent variable (i.e. density or occupancy) in the corresponding time interval. All
#                              values in this column must be non-negative, and there must be at least one value that is non-zero.
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
#
# Requirements:
#
#   R (Version >= 3.6.1)
#   R Package: data.table (Version >= 1.12.8)
#   R Package: gamlss (Version >= 5.1-5)
#   R Package: ggplot2 (Version >= 3.2.1)
#   R Package: ggpubr (Version >= 0.2.5)
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
if (nargs < 9) {
  cat('ERROR - Too few command-line arguments specified...\n')
  q(save = 'no', status = 1)
} else if (nargs > 9) {
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
}

# Perform checks on the command-line arguments and report their values
cat('\n')
cat('Path to R modules:                                  ', path_to_modules, '\n')
if (!dir.exists(path_to_modules)) {
  cat('ERROR - The directory path to the required R modules does not exist...\n')
  q(save = 'no', status = 1)
}
cat('Input data file:                                    ', data_file, '\n')
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
acceptable_values = c('FF', 'GS1935', 'GS1935kjf', 'GB1959', 'GB1959kjf', 'ED1961', 'ED1961kjf')               #### FINISH
if (!is.element(functional_form_model, acceptable_values)) {
  cat('ERROR - The command-line argument "functional_form_model" does not have an acceptable value...\n')
  q(save = 'no', status = 1)
}
cat('Functional form model:                              ', functional_form_model, '\n')
acceptable_values = c('GCV', 'gaussian_quadratic_var', 'gaussian_bspline_var')                                 #### FINISH
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
cat('Loading some R functions specific to this script...\n')
tryCatch(
  { source(file.path(path_to_modules, 'fitfun_functions.R')) },
  error = function(cond) { cat('ERROR - Failed to load the R functions...\n')
                           q(save = 'no', status = 1) }
)

# Define the names of the output files
output_file1 = file.path(output_dir, paste0('Fit.Summary.', fd_type, '.', functional_form_model, '.', noise_model, '.txt'))
output_file2 = file.path(output_dir, paste0('Fit.Curves.', fd_type, '.', functional_form_model, '.', noise_model, '.txt'))
output_file3 = file.path(output_dir, paste0('Fit.Predictions.', fd_type, '.', functional_form_model, '.', noise_model, '.txt'))
#output_file4 = ???                                                                                                           #### FINISH
#output_file5 = ???                                                                                                           #### FINISH

# If the output directory "output_dir" already exists
if (dir.exists(output_dir)) {

  # Report
  cat('\n')
  cat('The output directory already exists:', output_dir, '\n')
  cat('Checking if any of the output files already exist...\n')

  # If the output file "Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt" already exists
  if (file.exists(output_file1)) {

    # If the command-line argument "overwrite" is set to 'yes'
    if (overwrite == 'yes') {

      # Remove the output file "Fit.Summary.<fd_type>.<functional_form_model>.<noise_model>.txt"
      cat('Removing the output file:', output_file1, '\n')
      tryCatch(
        { file.remove(output_file1) },
        error = function(cond) { cat('ERROR - Failed to remove the output file...\n')
                                 q(save = 'no', status = 1) }
      )

    # If the command-line argument "overwrite" is set to 'no'
    } else {

      # Stop the script without doing anything
      cat('ERROR - The following output file already exists:', output_file1, '\n')
      q(save = 'no', status = 1)
    }
  }

  # If the output file "Fit.Curves.<fd_type>.<functional_form_model>.<noise_model>.txt" already exists
  if (file.exists(output_file2)) {

    # If the command-line argument "overwrite" is set to 'yes'
    if (overwrite == 'yes') {

      # Remove the output file "Fit.Curves.<fd_type>.<functional_form_model>.<noise_model>.txt"
      cat('Removing the output file:', output_file2, '\n')
      tryCatch(
        { file.remove(output_file2) },
        error = function(cond) { cat('ERROR - Failed to remove the output file...\n')
                                 q(save = 'no', status = 1) }
      )

    # If the command-line argument "overwrite" is set to 'no'
    } else {

      # Stop the script without doing anything
      cat('ERROR - The following output file already exists:', output_file2, '\n')
      q(save = 'no', status = 1)
    }
  }

  # If the output file "Fit.Predictions.<fd_type>.<functional_form_model>.<noise_model>.txt" already exists
  if (file.exists(output_file3)) {

    # If the command-line argument "overwrite" is set to 'yes'
    if (overwrite == 'yes') {

      # Remove the output file "Fit.Predictions.<fd_type>.<functional_form_model>.<noise_model>.txt"
      cat('Removing the output file:', output_file3, '\n')
      tryCatch(
        { file.remove(output_file3) },
        error = function(cond) { cat('ERROR - Failed to remove the output file...\n')
                                 q(save = 'no', status = 1) }
      )

    # If the command-line argument "overwrite" is set to 'no'
    } else {

      # Stop the script without doing anything
      cat('ERROR - The following output file already exists:', output_file3, '\n')
      q(save = 'no', status = 1)
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

# Check that the data file "data_file" exists
cat('\n')
cat('Reading in the data file:', data_file, '\n')
if (!file.exists(data_file)) {
  cat('ERROR - The data file does not exist...\n')
  q(save = 'no', status = 1)
}

# Read in the data file "data_file"
tryCatch(
  { data = fread(data_file, header = FALSE) },
  error = function(cond) { cat('ERROR - Failed to read in the data file...\n')
                           q(save = 'no', status = 1) }
)

# Perform basic checks on the data that have been read in
ndata = nrow(data)
if (ndata == 0) {
  cat('ERROR - The data file is empty...\n')
  q(save = 'no', status = 1)
}
if (ncol(data) != 3) {
  cat('ERROR - The data file does not have exactly three columns...\n')
  q(save = 'no', status = 1)
}
if (!is.numeric(data$V1)) {
  cat('ERROR - The data in the first column of the data file are not numeric...\n')
  q(save = 'no', status = 1)
}
if (!is.numeric(data$V2)) {
  cat('ERROR - The data in the second column of the data file are not numeric...\n')
  q(save = 'no', status = 1)
}
if (!is.numeric(data$V3)) {
  cat('ERROR - The data in the third column of the data file are not numeric...\n')
  q(save = 'no', status = 1)
}
if (ndata < 5) {
  cat('ERROR - The data file has less than 5 data lines...\n')
  q(save = 'no', status = 1)
}
if (!all(is.finite(data$V1))) {
  cat('ERROR - At least one data value in the first column is infinite...\n')
  q(save = 'no', status = 1)
}
if (!all(is.finite(data$V2))) {
  cat('ERROR - At least one data value in the second column is infinite...\n')
  q(save = 'no', status = 1)
}
if (!all(is.finite(data$V3))) {
  cat('ERROR - At least one data value in the third column is infinite...\n')
  q(save = 'no', status = 1)
}
if (any(data$V2 < 0.0)) {
  cat('ERROR - At least one data value in the second column is negative...\n')
  q(save = 'no', status = 1)
}
if (all(data$V2 == 0.0)) {
  cat('ERROR - All of the data values in the second column are zero...\n')
  q(save = 'no', status = 1)
}
if (any(data$V3 < 0.0)) {
  cat('ERROR - At least one data value in the third column is negative...\n')
  q(save = 'no', status = 1)
}
if (all(data$V3 == 0.0)) {
  cat('ERROR - All of the data values in the third column are zero...\n')
  q(save = 'no', status = 1)
}
cat('No. of lines read in:    ', ndata, '\n')

# Make a further check on "upper_density"
if (upper_density < max(data$V2)) {
  cat('ERROR - The command-line argument "upper_density" is less than the maximum value of the density data...\n')
  q(save = 'no', status = 1)
}


################################################################################################################################################
# Flow-density fundamental diagrams

# If the form of the FD relationship is flow-density
if (fd_type == 'Flow.Density') {

  # If the model component for the functional form of the flow-density relationship is the free-flow model (FF)
  if (functional_form_model == 'FF') {

    # If the model component for the noise in the flow-density relationship is defined as independent observations that follow a Gaussian distribution
    # with constant variance (GCV)
    if (noise_model == 'GCV') {

      # Load the required R module for performing the fit
      cat('\n')
      cat('Calling the required R module for performing the fit...\n')
      tryCatch(
        { source(file.path(path_to_modules, 'fit_flow_density_with_FF_GCV.R')) },
        error = function(cond) { cat('ERROR - Failed to load the required R module...\n')
                                 q(save = 'no', status = 1) }
      )

      # Fit the chosen GAMLSS model to the data                                      #### FINISH CLEANUPS IN THIS SECTION
      tryCatch(
        { model_obj = fit_flow_density_with_FF_GCV(data, ngrid, upper_density, output_file1, output_file2, output_file3) },
        error = function(cond) { cat('ERROR - Failed to fit the GAMLSS model for unknown reasons...\n')
                                 if (file.exists(output_file1)) { file.remove(output_file1) }
                                 if (file.exists(output_file2)) { file.remove(output_file2) }
                                 if (file.exists(output_file3)) { file.remove(output_file3) }
                                 q(save = 'no', status = 1) }
      )
    }
  }


################################################################################################################################################
# Speed-density fundamental diagrams

# If the form of the FD relationship is speed-density
}

# Finish
cat('\n')
cat('Finished!\n')
