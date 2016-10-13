# prophetic-granger-causality

Welcome to the Prophetic Granger Causality repo! The paper that this repo supports is:

Prophetic Granger Causality to infer Gene Regulatory Networks

Daniel E. Carlin, Evan O. Paull, Kiley Graim, Chris Wong, Adrian Bivol, Peter Ryabinin, Kyle Ellrott, Artem Sokolov, Joshua M. Stuart

If you use the GENIE3 code, please site:

Inferring regulatory networks from expression data using tree-based methods
Van Anh HUYNH-THU, Alexandre IRRTHUM, Louis WEHENKEL, Pierre GEURTS
PLoS ONE vol. 5(9): e12776

Prereqs:

You will have to install R (http://www.r-project.org).

For any of the GENIE3 variants, you must also install the R package randomForest.
Random Forests references:
Breiman L (2001) Random forests. Machine Learning 45: 5-32
http://stat-www.berkeley.edu/users/breiman/RandomForests/

This installation can be done from R with this command (as root under Linux):
> install.packages("randomForest")

You will nee python 2.7 for the prior heat kernel calculations, along with these packages:
numpy
scipy
optparse

Usage:

The complete dream8 winning submission can be recreated with:

make dream8

Within R, the data methods take a single data cube as input.  This cube Z is a list of n by m
matrices, where each element of the list is a replicate of n time steps by m genes.  The result that
is returned is a single m by m matrix of gene-gene interaction predictions.

source('dream8-granger-full.R')
source('prophetic-Genie3.R')

#for PGC

PGC_result=full.granger(Z)

#for Prophetic GENIE3

PG3_result=prophetic_GENIE3(Z)