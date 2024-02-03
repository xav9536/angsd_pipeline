

argv <- commandArgs(T)
DIST = as.numeric(argv[1])/2

library(data.table)
library(GenomicRanges)

fai = read.table("02_info/genome.fasta.fai")

region_num = as.character(read.table("02_info/regions_number.txt", colClasses = "character")[,1])

for(NUM in region_num){

LR_SNP = fread(paste0("03A_ngsparalog/all_maf0.01_pctind0.75_maxdepth4_chr", NUM, ".ngsparalog"))
LR_SNP$p.value = p.adjust(0.5*pchisq(LR_SNP$V5,df=1,lower.tail=FALSE), method="BH")
LR_SNP$deviant = LR_SNP$p < 0.001

dev = data.frame(mid = LR_SNP$V2[LR_SNP$deviant == TRUE])
dev$start = dev$mid - DIST
dev$end = dev$mid + DIST
dev$chr = LR_SNP$V1[1]
dev$start[dev$start < 1] = 1 # Correct for region "starting" before the beggining of the chromosome

gdev = makeGRangesFromDataFrame(dev, seqnames.field = "chr", start.field = "start", end.field = "end")
gdevu = union(gdev, gdev)
devu = as.data.frame(gdevu)
gdevus = gdevu[width(gdevu) > 151,] #Don't mask isolated deviant SNPs
devus = as.data.frame(gdevus)
devus[devus$end > fai[fai$V1 == devus$seqnames[1], "V2"], "end"] = fai[fai$V1 == devus$seqnames[1], "V2"] # Correct for region "ending" before the beggining of the chromosome

print(paste(round((sum(devus$width)/fai[fai$V1 == devus$seqnames[1], "V2"])*100, digits = 1), "% of chr", NUM, "masked"))

mask = devus

write.table(mask[,1:3],
            file = paste0("02_info/mask_by_chr/mask_deviant_chr", NUM, ".bed"),
            quote = F, row.names = F, col.names = F, sep = "\t")

}
