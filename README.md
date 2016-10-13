# prophetic-granger-causality

Welcome to the Prophetic Granger Causality repo! This code allows you to recreate the results of experiments presented in the work by Carlin, *et. al.* [1]

If you use the GENIE3 portion of the code, please cite the corresponding method [2].

## Prerequisites

You will have to install R (http://www.r-project.org).

For any of the GENIE3 variants, you must also install the R package randomForest [1]. This installation can be done from R with this command (as root under Linux):
```
> install.packages("randomForest")
```

You will need python 2.7 for the prior heat kernel calculations, along with these packages:
- numpy
- scipy
- optparse

## Usage

The complete dream8 winning submission can be recreated with:
```
make dream8
```
To use PGC, you will need to at least do:
```
make lasso.so
```
Within R, the data methods take a single data cube as input.  This cube `Z` is a list of `n` by `m`
matrices, where each element of the list is a replicate of `n` time steps by `m` genes.  The result that
is returned is a single `m` by `m` matrix of gene-gene interaction predictions.

Begin by executing the following:
```
source('dream8-granger-full.R')
source('prophetic-Genie3.R')
```
From here, if you are interested in PGC, use the following:

```
PGC_result=full.granger(Z)
```
To use Prophetic GENIE3, use the the following instead:
```
PG3_result=prophetic_GENIE3(Z)
```
See also the synapse for this project:

https://www.synapse.org/#!Synapse:syn2347433/wiki/62276

## Citations

[1] Daniel E. Carlin, Evan O. Paull, Kiley Graim, Chris Wong, Adrian Bivol, Peter Ryabinin, Kyle Ellrott, Artem Sokolov, Joshua M. Stuart. "Prophetic Granger Causality to infer Gene Regulatory Networks."

[2] Van Anh HUYNH-THU, Alexandre IRRTHUM, Louis WEHENKEL, Pierre GEURTS. "Inferring regulatory networks from expression data using tree-based methods." *PLoS ONE vol. 5(9): e12776*

[3] Breiman L (2001) "Random forests." *Machine Learning 45: 5-32.*
http://stat-www.berkeley.edu/users/breiman/RandomForests/
