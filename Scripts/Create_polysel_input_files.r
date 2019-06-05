library(data.table)
library(matrixStats)
library(foreach)

#set working directory
work.dir <- "~/Dropbox/Selection_test/"
setwd(work.dir)

#path to raw files
raw.path <- file.path(work.dir, "Raw_files")

#read in genes and go categories
genes <- fread(file.path(raw.path, "Entrez_Gene_IDs_positions_HG19.txt"))
go.cats <- fread(file.path(raw.path, "PolySel_GO_pathways_HG19.txt"))
go.cats[, setID:=as.numeric(factor(Pathway))]

#read in corrected Z scores
ll <- list.files(file.path(work.dir, "GeneScores"))
file.in <- ll[grep(".*GeneLengthCorrected.txt", ll)]
scores <- fread(file.path(work.dir, "GeneScores", file.in))

#output file path
out.path <- file.path(work.dir, "PolySel/Input")

pops <- unique(scores$Population)
for(pop in pops){
   print(pop)
   mm1 <- merge(scores[Population==pop, .(GeneID, objScore, Chr, Start, End)],
                go.cats, by="GeneID")
   mm2 <- merge(mm1, genes[, .(GeneID)], by="GeneID")

   mm2[, objID:=as.numeric(factor(GeneID))]
   obj.info <- mm2[!duplicated(GeneID),
                   .(objID, objStat=objScore, objName=GeneID, GeneLength=End-Start,
                     chr=Chr, startpos=Start, endpos=End, strand=NA)]

   fwrite(obj.info, file.path(out.path, paste0(pop, ".ObjInfo.txt")), sep="\t", na="NA", quote=F)

   set.info <- mm2[!duplicated(Pathway)][order(Source, Pathway),
                                         .(setID, setName=Pathway, setSource=Source)]
   fwrite(set.info, file.path(out.path, paste0(pop, ".SetInfo.txt")), sep="\t", na="NA", quote=F)

   set.obj <- mm2[order(setID, objID), .(setID, objID)]
   fwrite(set.obj, file.path(out.path, paste0(pop, ".SetObj.txt")), sep="\t", na="NA", quote=F)
}

