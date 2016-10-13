## dream8-prior.R
##
## by Artem Sokolov

pid2ker <- function()
  {
    fns <- list.files( "prior-pid" )
    for( fn in fns )
      {
        cat( "Processing", fn, "\n" )
        
        ## Convert a .pid to a .sif
        fn.pid <- paste( "prior-pid", fn, sep="/" )
        cmd <- paste0( "./pid2sif.py < ", fn.pid, " > temp.sif" )
        system( cmd )

        ## Apply the diffusion kernel to the .sif
        fn.ker <- paste( "prior-ker", fn, sep="/" )
        cmd <- paste0( "./makeKernel -n temp.sif -o ", fn.ker )
        system( cmd )

        ## Clean up
        cmd <- "rm temp.sif"
        system( cmd )
      }
  }

gen.prior <- function()
  {
    ## Retrieve the identifiers
    gid <- read.table( "global-idmap.tab", sep="\t" )
    hugo <- unique( as.character(gid[,1]) )
    n <- length( hugo )

    ## Initialize the results matrix
    prior <- matrix( 0, n, n )
    colnames(prior) <- hugo
    rownames(prior) <- hugo

    ## Load the kernel matrices
    fns <- paste( "prior-ker", list.files( "prior-ker" ), sep="/" )
    for( fn in fns )
      {
        ## Load the matrix
        X <- read.delim( fn, row.names=1, check.names=FALSE )
        X <- as.matrix(X)
        diag(X) <- 0
        stopifnot( all( rownames(X) == colnames(X) ) )

        ## Select the HUGO submatrix
        nn <- intersect( rownames(X), hugo )
        X <- X[nn,nn]

        if( sum(X) == 0 ) next

        prior[nn,nn] <- X
      }

    ## Normalize the evidence scores
    prior <- prior / max(prior)

    ## Cell and stimulus names
    cname <- c( "BT20", "BT549", "MCF7", "UACC812" )
    sname <- c( "EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1", "PBS", "Serum" )

    ## Traverse the cell lines
    clres <- list()
    for( clname in cname )
      {
        cat( "Processing", clname, "\n" )
        
        ## Retrieve the set of anti-body ids associated with this cell
        fnid <- paste( "data/", clname, "-idmap.tab", sep="" )
        idmap <- read.table( fnid, header=TRUE )
        nms <- as.character( idmap[,1] )
        n <- length(nms)

        res1 <- matrix( 0, n, n )
        rownames(res1) <- nms
        colnames(res1) <- nms
        
        for( rnm in nms )
          {
            ## Retrieve the row names
            rj <- which( gid[,2] == rnm )
            if( length(rj) == 0 ) next
            rh <- intersect( as.character(gid[rj,1]), rownames(prior) )
            if( length(rh) == 0 ) next
            
            for( cnm in nms )
              {
                cj <- which( gid[,2] == cnm )
                if( length(cj) == 0 ) next
                ch <- intersect( as.character(gid[cj,1]), colnames(prior) )
                if( length(ch) == 0 ) next
                
                res1[rnm, cnm] <- mean( prior[rh,ch] )
              }
          }

        clres[[clname]] <- res1 / max(res1)
      }

    res <- list()
    for( cn in cname )
      for( sn in sname )
        {
          nm <- paste( cn, sn, sep="." )
          res[[nm]] <- clres[[cn]]
        }

    res
  }

