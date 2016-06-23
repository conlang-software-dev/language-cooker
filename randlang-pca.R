library("FactoMineR")

load("bigwals-sample.RData")

data <- read.csv("bigwals-data-imputed-mca.csv", row.names = 1, check.names = F)
cltreedata <- as.matrix(read.csv("bigwals-cltree.csv",row.names = 1))
cltree <- graph_from_edgelist(cltreedata[,1:2])
E(cltree)$weight <- as.numeric(cltreedata[,3])

res <- PCA(data,T,Inf,graph = F)
vars <- apply(data,2,var)
evecs <- sweep(res$svd$V,1,vars,'*')

data <- data[result,]

randlang <- round(colMeans(data) + rowSums(sweep(evecs,2,sapply(res$eig$eigenvalue,function(x) rnorm(1,sd = sqrt(x))),"*")))
randlang <- as.matrix(pmin(pmax(randlang,0),apply(data,2,max)))

randlang

langdiff <- rowMeans(abs(t(randlang)/apply(data,2,max,na.rm=T) - (data/apply(data,2,max,na.rm=T)))*E(cltree)$weight, na.rm = T)
message(paste("Your language is closest to ",names(which.min(langdiff))," (distance = ",min(langdiff),").",sep = ""))