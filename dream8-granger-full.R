## dream8-granger-full.R - Granger-based predictor that utilizes all the data
##
## by Artem Sokolov

dyn.load( "lasso.so" )

## Solves the LASSO problem: argmin (y - Xw)^2 + \lambda |w|
optLASSO <- function( X, y, lambda, max.iter=1000, eps=1e-10 )
  {
    n <- nrow(X)
    p <- ncol(X)

    ## Set the initial parameter estimates
    w <- rep( 0, p )
    b <- mean( y )
    S <- X %*% w + b

    ## Call the C routine
    res <- .C( "optLASSO", as.double(X), as.double(y), 
              as.double(lambda), as.double(S), as.integer(n), as.integer(p),
              as.integer(max.iter), as.double(eps),
              w = as.double( w ), b = as.double(b), as.integer(1) )
    names( res$w ) <- colnames(X)
    
    list( w = res$w, b = res$b )
  }

## Solves the LASSO problem: argmin (y - Xw)^2 + \lambda |w|
## The meta-parameter \lambda is chosen such that w[ipvt] are barely zero
granger.lasso <- function( X, y, ipvt )
  {
    ## Compute the starting limits
    b <- mean(y)
    ll.max <- max( abs(apply( (y-b)*X, 2, mean )) )
    ll.min <- 0

    colnames(X) <- NULL
    
    ## Binary search for the proper \lambda parameter
    for( iter in 1:20 )
      {
        ll.cur <- mean( c(ll.max, ll.min) )
        m <- optLASSO( X, y, ll.cur )
        if( sum(m$w[ipvt]) == 0 ) ll.max <- ll.cur
        else ll.min <- ll.cur
      }

    ## Auto-regression may be slightly off zero depending on where binary
    ##     search terminated
    m$w[ipvt] <- 0
    m$w
  }

full.granger <- function( Zl )
  {
    stopifnot( length(Zl) > 1 )
    
    ## Verify consistency: number of proteins
    p <- unique(unlist( lapply( Zl, ncol ) ))
    stopifnot( length(p) == 1 )

    ## Verify consistency: the order of proteins
    vv <- lapply( Zl, function(z) { all( colnames(z) == colnames(Zl[[1]]) ) } )
    stopifnot( all(unlist(vv)) == TRUE )

    ## Reduce to a common set of timepoints
    pnms <- colnames( Zl[[1]] )
    rnms <- rownames( Zl[[1]] )
    invisible(lapply( Zl, function(z) { rnms <<- intersect(rnms, rownames(z)) } ))
    Zl <- lapply( Zl, function(z) {z[rnms,]} )
    stopifnot( length( rnms ) > 1 )

    ## Skip predicting the "baseline" point
    tstart <- 1
    if( rnms[1] == "bl" ) tstart <- 2

    ## Causal adjacency matrices
    A.past <- matrix( 0, p, p )
    rownames( A.past ) <- pnms
    colnames( A.past ) <- pnms
    A.fut <- A.past
    
    ## Traverse the targets
    for( jp in 1:p )		## Index j over the proteins
      {
        cat( "Processing", colnames(Zl[[1]])[jp], "\n" )
        yl <- lapply( Zl, function(z) {z[,jp]} )
        Xl <- lapply( Zl, function(z) {z[,-jp]} )

        ## Build the input space
        ##   use transpose of Xl[[i]]'s to keep timepoints together
        X.nonAR <- do.call( rbind, lapply( Xl, function(x1) {c(t(x1))} ) )
        n.nonAR <- ncol(X.nonAR)

        ## Keep track of (past) -> (present) and (present) -> (future) links
        W.past <- matrix( 0, 0, (p-1) )
        colnames( W.past ) <- colnames( Xl[[1]] )
        W.fut <- W.past
        
        ## Traverse the timepoints
        for( it in tstart:length(rnms) )	## Index i over the timepoints
          {
            cat( "Timepoint", it, "\n" )
            
            ## Build the response
            y <- unlist( lapply( yl, function(y1) {y1[it]} ) )

            ## Determine the index threshold, past which the dependence is
            ##  considered to be from the future, causing a reverse in the
            ##  causal link
            fut.thresh <- it * (p-1)

            ## Append the auto-regression component
            X.AR <- do.call( rbind, lapply( yl, function(y1) {y1[-it]} ) )
            X <- cbind( X.nonAR, X.AR )
            ipvt <- (n.nonAR+1):ncol(X)		## AR pivot

            ## Solve the regression problem using LASSO regularization
            w <- granger.lasso( X, y, ipvt )
            wsum <- sum( abs(w) )

            ## Normalize the weights to [-1, 1]
            if( wsum > 0 )
              w <- w / wsum

            ## Split the weights into past and future causality
            w.past <- w[1:fut.thresh]
            w.fut <- numeric(0)
            if( fut.thresh < n.nonAR )
              w.fut <- w[(fut.thresh+1):n.nonAR]

            ## Compute the average contribution of each time series
            W.past <- rbind( W.past, apply( matrix(w.past, nrow=p-1),1,mean ) )
            if( length( w.fut ) > 0 )
              W.fut <- rbind( W.fut, apply( matrix(w.fut, nrow=p-1),1,mean ) )
          }

        ## Entry (i, j) specifies a causal link from protein i to protein j
        ##  So, "past" entries go into the jth column, while "future" entries
        ##  go into the jth row
        A.past[colnames(W.past),jp] <- apply( W.past, 2, mean )
        A.fut[jp,colnames(W.fut)] <- apply( W.fut, 2, mean )
      }

    ## Do a simple average between past and future
    (A.past + A.fut) / 2
  }

combine.reps <- function( X, fcomb=median )
  {
    ## Retrieve the timepoint names
    nms <- unlist(lapply( strsplit( rownames(X), "-" ), "[[", 1 ))
    nml <- unique( nms )

    ## Combine the replicates
    res <- matrix( 0, length(nml), ncol(X) )
    rownames( res ) <- nml
    colnames( res ) <- colnames(X)
    for( nm in nml )
      res[nm,] <- apply( X[which( nms == nm ),,drop=FALSE], 2, fcomb )

    res
  }

process.cell <- function( clname )
  {
    ## Stimuli
    stims <- c( "EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1", "PBS", "Serum" )

    ## Load the data
    fn <- paste( "data/", clname, "-main.RData", sep="" )
    load( fn )

    ## Traverse the stimuli
    res <- list()
    for( s in stims )
      {
        cat( "Processing", clname, s, "\n" )
        nm <- paste( clname, s, sep="." )
        X <- lapply( get(nm), combine.reps )
        res[[nm]] <- full.granger( X )
      }

    res
  }

main.granger <- function()
  {
    res1 <- process.cell( "BT20" )
    res2 <- process.cell( "BT549" )
    res3 <- process.cell( "MCF7" )
    res4 <- process.cell( "UACC812" )

    res <- c( res1, res2, res3, res4 )
    res
  }

