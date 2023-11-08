argv <- commandArgs(T)
file<-argv[1]
OUTFILE_sites <- argv[2]
threshold<-as.numeric(argv[3])


ngsP = read.table(file)

ngsP$pval <- 0.5*pchisq(ngsP$V5,df=1,lower.tail=FALSE) # append column of p-values
ngsP$pval.adj <- p.adjust(ngsP$pval, method="BH")

ngsP_nonparalog = ngsP[ngsP$pval.adj > threshold,]
ngsP_paralog = ngsP[ngsP$pval.adj <= threshold,]

write.table(ngsP_nonparalog[,1:2],
              paste0(OUTFILE_sites, "_canonical"),
              quote = F, row.names = F, col.names = F)
write.table(ngsP_paralog[,1:2],
              paste0(OUTFILE_sites, "_deviant"),
              quote = F, row.names = F, col.names = F)
