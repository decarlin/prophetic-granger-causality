

source( "dream8-data.R" )
load.all()

dream8_test_dat<-read.table('validation_data/dream8_test_data.txt', header=FALSE)
colnames(dream8_test_dat)<-c('CellLines','Ligand','Probe','RespondedToMTorInhibition')

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

load.cell <- function( clname )
{
    ## Stimuli
    stims <- c( "EGF", "FGF1", "HGF", "IGF1", "Insulin", "NRG1", "PBS", "Serum" )

    ## Load the data
    fn <- paste( "data/", clname, "-main.RData", sep="" )
    load( fn )

    data_out<- list()

    for( s in stims )
      {
        cat( "Processing", clname, s, "\n" )
        nm <- paste( clname, s, sep="." )
        X <- lapply( get(nm), combine.reps )
	data_out[[nm]]<-X
	}

data_out

}

BT20_train_data<-load.cell("BT20")
BT549_train_data<-load.cell("BT549")
MCF7_train_data<-load.cell("MCF7")
UACC812_train_data<-load.cell("UACC812")

save(list = c('BT20_train_data','BT549_train_data','MCF7_train_data','UACC812_train_data','dream8_test_dat'), file = "dream8_test_and_train.RData")