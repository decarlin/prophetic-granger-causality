dream8 : lasso.so
	R --vanilla --slave < main.R

lasso.so : lasso.c
	R CMD SHLIB lasso.c
