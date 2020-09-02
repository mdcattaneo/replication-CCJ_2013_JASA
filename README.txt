######################################################################
## Empirical Illustration Replication Files.
##
## Cattaneo, Crump and Jansson (2013): "Generalized Jackknife 
## Estimators of Weighted Average Derivatives (with discussions and 
## rejoinder)", Journal of the American Statistical Association 
## 108(504): 1243-1268, December 2013.
##
## Contact email: cattaneo@umich.edu
##
######################################################################

MAIN DESCRIPTION:

The main file is "replic-empillus.R", while the other files are all
called by this file.

1) Put all files in the same directory, and change the setwd(.)
command in "replic-empillus.R" accordingly.

2) Note that I am giving you the windows object file of our C routines
"C_Kernel_WAD-jack.win64.so". This works in windows machines. (If you
are working on a Linux machine, then use the other object file
"C_Kernel_WAD-jack.linux64.so".)

3) Source C code is also included. A few specific comments:

3.a) This code is meant for continuous covariates, as assumed in the
paper.

3.b) The function "wad" estimates the first component of the weighted
average derivative and its standard error, for a given dataset and
bandwidth, using the Gaussian product kernel of the chosen order (2,
4, 6, 8 or 10). The code allows for an arbitrary number of covariates.

4) The main R code proceeds as follows:

Step 1) WAD is estimated for a grid of bandwidths around the selected
bandwidth (e.g., one of our ROT bandwidth choices, coded up in
"replic-empillus-functions.R" and used in "replic-empillus.R").

Step 2) The generalized jackknife estimator is constructed by using
the adjacent WAD estimates with the appropriate weights, as explained
in the paper.

For this second step, you need to use the R function genjack(.), which
is located in "replic-empillus-functions.R". This function requires
the special structure of the data supplied (i.e., WAD estimates, 
bandwidths used, etc.). This structure is available in the main R file
("replic-empillus.R"). Note that this function is quite flexible: it 
allows you to jackknife using any combination of "lower" estimates or 
"higher" estimates, that is, any WAD estimates for bandwidths smaller
or bandwidths larger than the bandwidth originally chosen.

The rest of the code is, hopefully, self-explanatory.

