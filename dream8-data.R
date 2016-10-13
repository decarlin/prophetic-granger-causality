## dream8.R - DREAM 8 challenge
##
## by Artem Sokolov

parse.ligand <- function( X, lig )
  {
    j <- which( X[,"Stimulus"] == lig )
    X1 <- X[j,]

    res <- list()
    for( inh in unique( X1[,"Inhibitor"] ) )
      {
        j1 <- which( X1[,"Inhibitor"] == inh )
        inh1 <- inh; if( inh1 == "" ) inh1 <- "None"

        rn <- as.character(X1[j1,"Timepoint"])
        rn <- paste( rn, rep( 1:3, length=length(rn) ), sep="-" )
        
        res[[inh1]] <- as.matrix( X1[j1,-(1:4)] )
        rownames( res[[inh1]] ) <- rn
      }

    res
  }

parse.pairs <- function( X )
  {
    ligs <- c( "loLIG1+loLIG2", "loLIG1+hiLIG2", "hiLIG1+loLIG2", "hiLIG1+hiLIG2" )

    res <- list()
    for( lig in ligs )
      {
        j <- which( X[,"Stimulus"] == lig )
        
        rn <- as.character( X[j,"Timepoint"] )
        rn <- paste( rn, rep( 1:3, length=length(rn) ), sep="-" )

        res[[lig]] <- as.matrix( X[j,-(1:4)] )
        rownames( res[[lig]] ) <- rn
      }

    res
  }

load.insilico <- function()
  {
    fn <- "data/insilico.csv"
    X <- read.csv( fn, check.names=FALSE )

    InS.loLIG1 <- parse.ligand( X, "loLIG1" )
    InS.hiLIG1 <- parse.ligand( X, "hiLIG1" )
    InS.loLIG2 <- parse.ligand( X, "loLIG2" )
    InS.hiLIG2 <- parse.ligand( X, "hiLIG2" )

    InS.pairs <- parse.pairs( X )

    save( InS.loLIG1, InS.hiLIG1, InS.loLIG2, InS.hiLIG2, InS.pairs,
         file = "data/insilico.RData" )
  }

load.raw <- function( fn )
  {
    ## Load raw data
    X <- read.csv( fn, check.names=FALSE, header=FALSE )
    p <- ncol(X)

    ## Count the number of header rows
    nhr <- min(which( as.character( X[,1] ) != "" ))

    ## Retrieve the first part of the header
    h1 <- t(X[nhr,])
    j <- min(which( h1 == "" ))
    h1 <- h1[1:(j-1)]
    nhc <- length(h1)	## Number of header columns

    ## Retrieve the antibody names
    abj <- which( X[,nhc] == "Antibody Name" )
    hgj <- which( X[,nhc] == "HUGO ID" )
    h2 <- c(t(X[abj,(nhc+1):p]))
    h3 <- c(t(X[hgj,(nhc+1):p]))
    names( h3 ) <- h2

    ## Retrieve the data
    X <- X[-(1:nhr),]
    colnames(X) <- c( h1, h2 )
    
    list( X=X, idmap=h3 )
  }

parse.baseline <- function( X )
  {
    ## Find the baseline entries
    j1 <- which( X[,"Stimulus"] == "" )
    j2 <- which( X[,"Timepoint"] == "0min" )
    stopifnot( all( j1 == j2 ) )

    ## Traverse the inhibitors
    X1 <- X[j1,]
    res <- list()
    for( ii in unique( X1[,"Inhibitor"] ) )
      {
        j <- which( X1[,"Inhibitor"] == ii )
        res[[ii]] <- as.matrix(X1[j,-(1:4)])
        storage.mode( res[[ii]] ) <- "numeric"
        rownames( res[[ii]] ) <- paste( "Replicate", 1:nrow(res[[ii]]) )
      }
    res
  }

handle.dupes <- function( X )
  {
    rn <- rownames(X)
    j <- c()
    v <- c()
    for( ii in c( "5min", "15min", "30min", "60min", "2hr", "4hr" ) )
      {
        j1 <- which( rn == ii )
        n1 <- length(j1)
        v <- c( v, paste( ii, 1:n1, sep="-" ) )
        j <- c( j, j1 )
      }

    res <- X[j,]
    rownames(res) <- v
    res
  }

parse.stimulus <- function( X, ss )
  {
    ## Find the stimulus entries
    js <- which( X[,"Stimulus"] == ss )
    X1 <- X[js,]

    ## Traverse the inhibitors
    res <- list()
    for( ii in unique( X1[,"Inhibitor"] ) )
      {
        j <- which( X1[,"Inhibitor"] == ii )
        rn <- as.character(X1[j,"Timepoint"])
        res[[ii]] <- as.matrix( X1[j,-(1:4)] )
        storage.mode( res[[ii]] ) <- "numeric"
        rownames( res[[ii]] ) <- rn

        if( length( unique( rn ) ) < length(rn) )
          {
            warning( "Handling duplicate timepoints for ", ss, "/", ii )
            res[[ii]] <- handle.dupes( res[[ii]] )
          }
      }
    res
  }

save.idmap <- function( idmap, fn )
  {
    ab <- names( idmap )
    hg <- idmap
    names( hg ) <- NULL
    idm <- data.frame( Antibody = ab, "HUGO ID" = hg )
    write.table( idm, file=fn, quote=FALSE, row.names=FALSE, sep="\t" )
  }

load.main <- function( cln )
  {
    fn <- paste( "data/", cln, "_main.csv", sep="" )
    dat <- load.raw( fn )
    X <- dat$X
    idmap <- dat$idmap
    save.idmap( idmap, paste( "data/", cln, "-idmap.tab", sep="" ) )

    ## Compose the baseline
    X.bl <- parse.baseline( X )
    v.bl <- paste( cln, ".baseline", sep="" )
    assign( v.bl, X.bl )
    
    ## Retrieve the list of stimuli
    stims <- as.character(unique( X[,"Stimulus"] ))
    j <- which( stims == "" )
    if( length(j) > 0 ) stims <- stims[-j]
    v.ss <- paste( cln, stims, sep="." )

    ## Traverse the stimuli
    for( ss in stims )
      {
        X.ss <- parse.stimulus( X, ss )
        v <- paste( cln, ss, sep="." )
        assign( v, X.ss )
      }

    save( list=c( v.bl, v.ss ), file=paste( "data/", cln, "-main.RData", sep="" ) )
  }

load.all <- function()
  {
    load.main( "BT20" )
    load.main( "BT549" )
    load.main( "MCF7" )
    load.main( "UACC812" )
  }
