## dream8-res.R - Functions to write out DREAM8 results
##
## by Artem Sokolov

exp.names <- function()
  {
    cname <- c( "BT20", "BT549", "MCF7", "UACC812" )
    sname <- c( "EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1", "PBS", "Serum" )
    a <- expand.grid( cname, sname )
    paste( a[,1], a[,2], sep="." )
  }

antibody.names <- function( cellname )
  {
    fn <- paste( "data/", cellname, "-main.RData", sep="" )
    load( fn )
    varname <- paste( cellname, "baseline", sep="." )
    colnames( get( varname )[["DMSO"]] )
  }

verify.result <- function( res )
  {
    ## Verify the existence of experiments
    nms <- exp.names()
    msng <- setdiff( nms, names(res) )
    if( length( msng ) > 0 )
      stop( "Missing the following entries:\n", paste( msng, collapse="\n" ) )

    ## Verify the existence of In Silico
    jInS <- grep( "Insilico", names(res) )
    if( length(jInS) == 0 )
      warning( "Missing In Silico" )

    ## Verify dimensionality
    jBT20 <- grep( "BT20", names(res) )
    if( all(unlist(lapply( res[jBT20], dim )) == 48) == FALSE )
      stop( "All of the BT20 entries must be 48x48" )
    
    jBT549 <- grep( "BT549", names(res) )
    if( all(unlist(lapply( res[jBT549], dim )) == 45) == FALSE )
      stop( "All of the BT549 entries must be 45x45" )
    
    jMCF7 <- grep( "MCF7", names(res) )
    if( all(unlist(lapply( res[jMCF7], dim )) == 41) == FALSE )
      stop( "All of the MCF7 entries must be 41x41" )
    
    jUACC812 <- grep( "UACC812", names(res) )
    if( all(unlist(lapply( res[jUACC812], dim )) == 45) == FALSE )
      stop( "All of the UACC812 entries must be 45x45" )
    
    ## Verify individual entries
    lapply( res, function( z )
           {
             if( is.null( rownames(z) ) )
               stop( "One of the results entries is missing rownames" )
             if( is.null( colnames(z) ) )
               stop( "One of the results entries is missing colnames" )
             if( all( rownames(z) == colnames(z) ) == FALSE )
               stop( "rownames and colnames must match" )
           } )

    invisible( NULL )
  }

write.sif.eda <- function( X, pfx )
  {
    ## Remove loops and normalize the data
    diag(X) <- 0
    X <- abs(X)
    stopifnot( max(X) > 0 )
    X <- X / max(X)

    Z <- cbind( expand.grid( rownames(X), colnames(X) ), c(X) )
    j <- which( Z[,3] == 0 )
    Z <- Z[-j,]

    ## Write the .sif
    sif1 <- paste( Z[,1], sign(Z[,3]), Z[,2], sep="\t", collapse="\n" )
    j <- which( rowSums( X ) + colSums( X ) == 0 )
    if( length(j) > 0 )
      {
        sif2 <- paste( names(j), collapse="\n" )
        sif1 <- paste( sif1, sif2, sep="\n" )
      }
    fnSif <- paste( pfx, "sif", sep="." )
    write( sif1, fnSif )

    ## Write the .eda
    eda0 <- "EdgeScore"
    eda1 <- paste( Z[,1], " (", sign(Z[,3]), ") ", Z[,2], " = ",
                  round( Z[,3], 5 ), sep="", collapse="\n" )
    fnEda <- paste( pfx, "eda", sep="." )
    write( paste( eda0, eda1, sep="\n" ), fnEda )
  }

load.eda <- function( fn, abnames )
  {
    X <- read.table( fn )
    n <- length( abnames )
    res <- matrix( 0, n, n )
    rownames(res) <- abnames
    colnames(res) <- abnames

    ## Traverse the anti-bodies
    for( ab in abnames )
      {
        j <- which( X[,1] == ab )
        if( length(j) == 0 ) next
        jab <- as.character( X[j,3] )
        res[ ab, jab ] <- X[j,5]
      }

    res
  }

eda2res <- function( pfx, fnOut )
  {
    cnms <- c( "BT20", "BT549", "MCF7", "UACC812" )
    snms <- c( "EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1", "PBS", "Serum" )

    res <- list()
    for( cnm in cnms )
      {
        ## Load the anti-body names associated with this cell
        abnames <- antibody.names( cnm )

        ## Traverse the stimuli
        for( snm in snms )
          {
            fn <- paste0( pfx, cnm, ".", snm, ".eda" )
            rnm <- paste( cnm, snm, sep="." )
            res[[rnm]] <- load.eda( fn, abnames )
          }
      }

    verify.result( res )
    save( res, file=fnOut )
  }

main.res <- function( fnRes = "pred/final.RData", pfxTeam = "DC_TDC" )
  {
    ## Load the results
    load( fnRes )
    if( exists( "res" ) == FALSE )
      stop( "No object 'res' in file ", fnRes )
    verify.result( res )
    
    ## Create a temporary subdirectory
    dnCur <- getwd()
    dnSave <- paste( pfxTeam, "temp", sep="-" )
    dir.create( file.path( dnCur, dnSave ), showWarnings = FALSE )
    setwd( file.path( dnCur, dnSave ) )

    ## Traverse individual 'res' entries
    nms <- exp.names()
    for( nm in nms )
      {
        pfx <- paste( pfxTeam, sub( "\\.", "-", nm ), "Network", sep="-" )
        write.sif.eda( res[[nm]], pfx )
      }

    ## Generate a 'writeup'
    fnWup <- paste( pfxTeam, "Network-Writeup.txt", sep="-" )
    cmd <- paste( "touch", fnWup )
    system( cmd )

    ## Zip up the current entries
    fnZip <- paste( pfxTeam, "Network.zip", sep="-" )
    cmd <- paste( "zip -r", fnZip, "*" )
    system( cmd )

    ## Write out in-silico, as necessary
    j <- which( names(res) == "Insilico" )
    if( length(j) > 0 )
      {
        pfx <- paste( pfxTeam, "Network", "Insilico", sep="-" )
        write.sif.eda( res[["Insilico"]], pfx )

        ## Generate a 'writeup' for in-silico
        fnWup <- paste( pfxTeam, "Network-Insilico-Writeup.txt", sep="-" )
        cmd <- paste( "touch", fnWup )
        system( cmd )

        ## Zip up the in-silico predictions
        fnZip <- paste( pfxTeam, "Network-Insilico.zip", sep="-" )
        cmd <- paste( "zip -r", fnZip, "*Insilico*" )
        system( cmd )
      }

    ## Move the zip files to the original directory
    cmd <- paste( "mv *.zip", dnCur )
    system( cmd )

    ## Restore the working directory
    setwd( dnCur )

    ## Cleanup
    cmd <- paste( "rm -r", file.path( dnCur, dnSave ) )
    system( cmd )
  }

