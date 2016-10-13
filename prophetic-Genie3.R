source('GENIE3_R/genie3.R')

prophetic_GENIE3<-function(cube)
{
	#cube is a list with elements of replicates, each of which is rows=time,cols=genes

	all.genes=colnames(cube[[1]])
	times=rownames(cube[[1]])
	out=matrix(0,nrow=length(all.genes),ncol=length(all.genes))
	rownames(out)=all.genes
	colnames(out)=all.genes
	in.mat = matrix(NA,nrow=length(all.genes)*length(times), ncol=length(cube))	

	rownames(in.mat)=c(outer(all.genes,times,function(x,y){paste(x,y,sep='_')}))
	row.genes=rep(all.genes,length(times))
	row.times=c()

	for (t in times)
	{
		row.times=c(row.times,rep(t,length(all.genes)))
	}

	#arrange in.mat such that we have (probe,time) in rows, replicates in columns

	for (rep in 1:length(cube))
	{
		for (t in times)
		{
			for (g in all.genes)
			{
				cur.name=paste(g,t,sep='_')
				in.mat[cur.name,rep]=cube[[rep]][t,g]
			}
		}
	}

	#here is the Genie3 call
	w.m <- get.weight.matrix(in.mat)

	names(row.times)<-rownames(w.m)
	names(row.genes)<-rownames(w.m)

	#now consolidate times and flip causality of future points.  flipCausality tells which ones

	flipCausality<-outer(row.times,row.times,function(x,y){x<y})

	for (source in all.genes)
	{
		for (target in all.genes)
		{
			if (source==target){next}
			geneMask<-outer(c(row.genes==source),c(row.genes==target))

			out[source,target]=sum(out[source,target],w.m[which(geneMask & !flipCausality, arr.ind=TRUE)])
			out[target,source]=sum(out[target,source],w.m[which(geneMask & flipCausality, arr.ind=TRUE)])
		}
	}
out
}