#!/usr/bin/env Rscript

##Example command line
#Rscript --verbose --vanilla ~/Dropbox/Selection_test/Scripts/polysel_CMD.r ~/Dropbox/Selection_test/Scripts/ PolySel/Input PolySel/Output Anatolia_EF 1000000 3

library("qvalue")
library("Matrix")
library("igraph")
library("data.table")
library("ggplot2")
library("foreach")
library("qvalue")
library("doParallel")

args <- commandArgs(trailingOnly=TRUE)

##test example to run from within R
#args <- c("~/Dropbox/Selection_test/", "PolySel/Input", "PolySel/Output", "EHG", 10000, 3)


##pathway to files
pth <- args[1]

##pathway to data
data.path <- file.path(pth, args[2])

##pathway to output
out.path <- file.path(pth, args[3])

##population for analysis
groupname <- args[4]

##number of permutations to calculate
emp.nruns <- as.numeric(args[5])

##number of cores to use for permutations
n.cores <- as.numeric(args[6])

##read in polysel R functions
code.path <- file.path(pth, "Scripts")
source(file.path(code.path,'polysel.R'))

##create folders for output
dir.create(paste0(out.path, "/", groupname))
results.path <- file.path(out.path, groupname, "results")
dir.create(results.path)


##-----------------------------------------##
##read in data
##-----------------------------------------##

print(paste("Inputting data files..."))
#minsetsize <- 10
raw <- ReadSetObjTables(in.path=data.path,
                        set.info.file=paste0(groupname, ".SetInfo.txt"),
                        set.obj.file=paste0(groupname, ".SetObj.txt"),
                        obj.info.file=paste0(groupname, ".ObjInfo.txt"),
                        minsetsize=1,
                        obj.in.set=F,
                        merge.similar.sets=T)

set.info <- raw$set.info
obj.info <- raw$obj.info
set.obj <- raw$set.obj
set.info.lnk <- raw$set.info.lnk

write.table(set.info.lnk, quote=FALSE, sep="\t", row.names=FALSE,
            file=file.path(results.path, "merged.pathways.txt"))


#remove na's
no.scores <- obj.info[is.na(obj.info$objStat),]$objID
obj.info <- obj.info[!(obj.info$objID %in% no.scores), ]
set.obj <- set.obj[!(set.obj$objID %in% no.scores), ]
obj.info <- obj.info[!is.na(obj.info$objStat),]


##-----------------------------------------##
##Run new enrichment test
##-----------------------------------------##

print(paste("Running enrichment test..."))

#permutation function
permute.data <- function(obj.info, n.chr, n.genes, gene.pos, chr.ord.now){
   # 1. string together chromosomes
   chr.ord.now <- sample(1:n.chr, n.chr)
   new.ord <- unlist(gene.pos[chr.ord.now])
   r1 <- obj.info[new.ord, objStat]

   # 2. rotate scores
   rotate.now <- sample(2:n.genes, 1)
   r1[c(rotate.now:n.genes, 1:(rotate.now-1))]
}

#sumstat calculation function
sum.stat <- function(set.obj, ID, perm){
   # 3. compute new sumstat scores
   mm.e <- merge(data.table(objID=ID, objStat=perm), set.obj, by="objID")[order(setID)]
   mm.e[, lapply(.SD, sum), .SDcols=grep("objStat", names(mm.e)), by="setID"]
}

#p.value calculation function
compute.p.val <- function(obs, exp){
   # 4. compute p values
   e <- exp[, setID:=NULL]
   rowSums(obs <= e)
}

#save compressed null distributions
null.bins <- function(exp, lower, upper, intervals){
   #5 save a smaller version of the null (rather than all empirical vales)
   cc <- c(-Inf, seq(lower, upper, intervals), Inf)
   t(apply(m.exp[, grep("objStat", names(m.exp)), with=F], 1,
           function(x) table(cut(x, cc))))
}

#timer function
timer <- function(time.point, emp.nruns, I){
   if((I-1000) %% time.point == 0 & I>1000){ # estimate time remaining
      d.now <- as.numeric(difftime(Sys.time(), s0, units="secs"))
      est.time <- (emp.nruns-I+1000)*(d.now/(I-1000))
      print(paste("Completed", I-1000, "iterations"))
      mins <- est.time %/% 60
      secs <- round(est.time %% 60, 0)
      if(mins>=60){
         hours <- mins %/% 60
         mins <- mins %% 60
         print(paste("Estimated time remaining:", hours, "h", mins, "m", secs, "s"))
      }else{
         if(mins==0){
            print(paste("Estimated time remaining:", secs, "s"))
         }else{
            print(paste("Estimated time remaining:", mins, "m", secs, "s"))
         }
      }
   }
}

#break up emp.nruns into specific size iteration chunks
#(improves run time; 1000 seems optimal for ~10k genes)
get.blocks <- function(emp.nruns, block.size){
   if(emp.nruns<block.size){
      rb <- emp.nruns
   }else{
      onek.blocks <- rep(block.size, emp.nruns %/% block.size)
      remainder <- emp.nruns %% block.size
      if(remainder>0){
         rb <- c(onek.blocks, remainder)
      }else{
         rb <- onek.blocks
      }
   }
   return(rb)
}

#convert to datatable format
setDT(set.info)
setDT(obj.info)
setDT(set.obj)

#set up parallel backend for foreach %dopar%
registerDoParallel(cores=n.cores)

# order obj.info table (needed to compare rotated values against)
obj.info <- obj.info[order(chr, startpos)]
ID <- obj.info$objID

write.table(obj.info, quote=FALSE, sep="\t", row.names=FALSE,
            file=file.path(results.path, "obj.info.full.txt"))

#get relevant variables
n.genes <- obj.info[, .N]
n.paths <- set.info[, .N]
n.chr <- obj.info[, length(unique(chr))]

#determine matrix position of each gene for each chromosome
gene.pos <- foreach(i=1:n.chr) %do% obj.info[, which(chr==i)]

# compute observed sumstat scores
mm.o <- merge(set.obj, obj.info[, .(objID, objStat)], by="objID")
mm.n <- mm.o[order(setID), .(N=length(unique(objID))),
             by=c("setID")]
m.obs <- mm.o[order(setID), .(SumStat=sum(objStat, na.rm=T)),
              by=c("setID")]

#housekeeping
rm(raw, set.info.lnk, mm.o)
gc()

# compute expected sumstat scores
run.blocks <- get.blocks(emp.nruns, 1000)
I=0
sig.tests <- rep(0, n.paths)
s0 <- Sys.time()
for(l in run.blocks){
   I=I+l
   timer(1000, emp.nruns, I)
   pp <- foreach(i=1:l) %dopar% {
      permute.data(obj.info, n.chr, n.genes, gene.pos, chr.ord.now)
   }
   #perm <- matrix(unlist(pp), ncol=l, byrow=TRUE)
   perm <- matrix(unlist(pp), ncol=l, byrow=FALSE)
   m.exp <- sum.stat(set.obj, ID, perm)
   sig.tests <- sig.tests + compute.p.val(obs=m.obs$SumStat, exp=m.exp)
   #sig.tests.fict <- sig.tests.fict + compute.p.val(obs=m.obs$SumStat.fict, exp=m.exp$Fictive)
   rm(perm, m.exp)
   gc()
}


##-----------------------------------------##
## Compute p and q values
##-----------------------------------------##

p.vals <- sig.tests/emp.nruns
q.vals <- qvalue(p.vals, pi0.method="smoother")

q.list <- list(P=p.vals, Q=q.vals)
save(q.list, file=file.path(results.path, "q.values.Rdata"))

##-----------------------------------------##
##write output (CONVERT TO POLYSEL FORMAT)
##-----------------------------------------##

out.file <- data.table(set.info[, .(setID, setSize=setSizeOrg)],
                       setScore=m.obs$SumStat, N=mm.n$N,
                       setP=p.vals, setQ=q.vals$qvalues,
                       set.info[, .(setName, setID.orig)])[order(setP)]

fwrite(out.file, quote=FALSE, sep="\t", row.names=FALSE,
       file=file.path(results.path, "setscores.txt"))



