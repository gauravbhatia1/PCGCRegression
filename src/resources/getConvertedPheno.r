### clean data
clean.data <- function(phen.file,fe.file){
	phen <- read.csv(phen.file,header=F,sep='',as.is=T)
	fe <- read.csv(fe.file,header=F,sep='',as.is=T)

	all <- union(fe[,2],phen[,2])
	good <- intersect(fe[,2],phen[,2])
	name.list = c();
	for(name in all){
		if(name %in% good){
			name.list <- c(name.list,name)
		}
	}

	# remove bad names from lists 
	fe <- fe[which(fe[,2] %in% name.list),]
	phen <- phen[which(phen[,2] %in% name.list),]
	

	#order the lists by phenotype order
        fe <- fe[match(phen[,2],fe[,2]),];
	if(!identical(fe[,2],phen[,2]))
	{
		print(head(fe[,2]));
		print(head(phen[,2]));
		print("must be identical");
		q();
	}
	
	ids = phen[,1:2];
	phen <- phen[,3]
	if(max(phen)==2){
		phen <- phen - 1
	}
	fe = data.matrix(fe[,-(1:2)]);
	if( fe.file == phen.file )
	{
		fe=rep(1,length(phen))
	}
	ret <- list(ids=ids,phen=phen,fe=fe)
	return(ret)
}

transformPheno <- function(phen,fe,K){
	# fit logistic regression 
	model <- glm(phen ~ fe , family="binomial")
	asc.probs <- model$fitted
	if(length(asc.probs) != length(phen)){
		# we have NAs 
		new.asc.probs <- array(NA,length(phen))
		new.asc.probs[which(!is.na(phen))] <- asc.probs
		asc.probs <- new.asc.probs 
	}
	P <- mean(phen,na.rm=T)
	asc.const <- (1-P)/P * K/(1-K)
	non.asc.probs <- asc.const * asc.probs / (1+asc.const * asc.probs - asc.probs)
	# convert probs to thresholds 
	non.asc.thresholds <- qnorm(1-non.asc.probs)
	# correct ascertainment
	phis <- dnorm(non.asc.thresholds)
	transpheno = (phen-asc.probs)/(sqrt(asc.probs*(1-asc.probs)));
	
	const1 <- K/(1-K) * (1-P)/P
	const2 <- (P-K)/(P*(1-K)) 

	x <- (1 - (asc.probs) * const2 )
	x <- x / sqrt(asc.probs*(1- asc.probs))
	x <- x * phis 
	x = x / (non.asc.probs*(1-const1)+const1)	
	transpheno = transpheno;

        # compute the population-wise variance of the thresholds using the law of total-variance
	non.asc.thresholds[non.asc.thresholds==Inf] = NA;
        var.cases <- var(non.asc.thresholds[phen==1],na.rm=T)
        var.controls <- var(non.asc.thresholds[phen==0],na.rm=T)
        mean.cases <- mean(non.asc.thresholds[phen==1],na.rm=T)
        mean.controls <- mean(non.asc.thresholds[phen==0],na.rm=T)
        total.var <- K*var.cases + (1-K)*var.controls + K*(1-K)*(mean.cases - mean.controls)^2

	ret <- list(phen=transpheno,multiplier=x,totalvar=total.var)
	return(ret);

}

##!/usr/local/bin/Rscript --vanilla
args <- commandArgs(trailingOnly = TRUE)

if(length(args)<5){
	stop("Usage: getConvertedPheno.R phen_file K fe_file transphenofile multiplierfile ")
}

# open and clean data
K <- as.double(args[2])
cd = clean.data(args[1],args[3]);
phenos_transformed = transformPheno(phen=cd$phen,fe=cd$fe,K=K);
outphenofilename <- args[4]
outmultfilename <- args[5]
transformed_pheno_data = cbind(cd$ids, phenos_transformed$phen);
write.table(transformed_pheno_data, file=outphenofilename, quote=F, row.names=F,col.names=F);
transformed_pheno_data = cbind(cd$ids, phenos_transformed$multiplier);
write.table(transformed_pheno_data, file=outmultfilename, quote=F, row.names=F,col.names=F);
cat(phenos_transformed$totalvar)
