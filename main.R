## main.R - Entry point of the code
##
## by Artem Sokolov

source( "dream8-data.R" )
source( "dream8-prior.R" )
source( "dream8-granger-full.R" )
source( "dream8-res.R" )
source( "dream8-nsmbl.R" )

main <- function()
  {
    ## Parse the data
    load.all()

    ## Run the prior
    pid2ker()
    res <- gen.prior()
    save( res, file="pred/prior.RData" )

    ## Run Granger Causality
    res <- main.granger()
    save( res, file="pred/granger.RData" )

    ## Generate the ensemble
    main.ensemble()

    ## Put the final predictions into the correct format
    main.res()
  }

main()

