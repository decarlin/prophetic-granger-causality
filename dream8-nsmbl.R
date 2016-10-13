## dream8-nsmbl.R - ensemble generation for the DREAM8 challenge
##
## by Artem Sokolov

loadRes <- function( fn )
  {
    load( fn )
    res
  }

make.ensemble <- function( fns, fcomb=mean, bInS=FALSE )
  {
    ## Load all the results files
    X <- list()
    for( fn in fns )
      X[[fn]] <- loadRes( fn )
    stopifnot( length(X) > 0 )

    ## The set of names to traverse
    nms <- exp.names()
    if( bInS == TRUE )
      nms <- c( nms, "Insilico" )
    
    ## Traverse experimental data
    res <- list()
    for( nm in nms )
      {
        ## Extract an normalize the entries
        X1 <- lapply( X, "[[", nm )
        X1 <- lapply( X1, abs )
        X1 <- lapply( X1, function(z) {z / max(z)} )

        ## Verify consistency
        nn <- unique(unlist(lapply( X1, length )))
        stopifnot( length(nn) == 1 )
        rn <- unique(unlist(lapply( X1, rownames )))
        rndiff <- unlist(lapply( X1, function(x1) {setdiff( rn, rownames(x1) )} ))
        stopifnot( length(rndiff) == 0 )

        ## Combine the entries
        n <- nrow(X1[[1]])
        A1 <- array( data=unlist(X1), dim = c(n,n,length(X1)) )
        A <- apply( A1, c(1,2), fcomb )

        ## Fix the row-/colnames
        rownames(A) <- rownames(X1[[1]])
        colnames(A) <- colnames(X1[[1]])
        res[[nm]] <- A
      }
    res
  }

main.ensemble <- function()
  {
    fns <- c( "pred/prior.RData",
             "pred/granger.RData" )

    res <- make.ensemble( fns )

    save( res, file="pred/final.RData" )
  }

