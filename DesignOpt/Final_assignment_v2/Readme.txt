# Project folder file for the course ME46060.
#
# The project is structured as follows:
# Main_optimizer_reduced is the file to be run for the optimization
# with reduced amount of variables. This file also uses a different
# finite difference optimizer. Supporting functions are (used for both
# optimizations:
#
# - GetInertia
# - Getshape(for plotting airfoils)
# - GetThickness
# - GetVolume
# - GetWeight
# - constraints
# - SensitivityAnalysis
#
# Both optimizers run the objective mask file to convert line search to
# an actual point in the design space for which analysis is carried out.
# This function then calls the real objective function. All other parameters
# are stored in opt_params.
#
# the file "ProblemAnalysis" was used to carry out a design space analysis
# and creating contour plots