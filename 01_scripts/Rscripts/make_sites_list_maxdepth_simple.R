#this R script uses the mafs file provided by the analyses on all individuals, with filter for quality and maf 0.05
#it simply extract the first columns with chr, position of each SNP and Major/minor alleles as determined in step 03
#output is a "sites" files that will allwo restraining subsequent analyses 
#(for instance maf by pop, FSt by pop pairs, etc to a limited number of loci

argv <- commandArgs(T)
INFILE <- argv[1]
OUTFILE_sites <- argv[2]

maf<-read.table(INFILE, header=T)
head(maf)
sites<-maf[,1:4]
sites_order <- sites[order(sites$chromo),] 

write.table(sites_order, OUTFILE_sites, row.names=F, col.names=F, sep="\t", quote=F)

