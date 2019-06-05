library(data.table)
library(ggplot2)
library(BEDMatrix)
library(parallel)

# function for making two pseudohaploid columns from single diploid vector


diploid2Pseudohaploid = function(x, name) {
  ID = !is.na(x)
  sorterColumn = rep(-9, length=length(x))
  sorterColumn[x==1 & ID] = rbinom(n = sum(x==1 & ID), size = 1, prob = 0.5) + 1
  firstColumn = ifelse(x==0 | x==2, x, ifelse(sorterColumn == 1, 0, 2))
  secondColumn = ifelse(x==0 | x==2, x, ifelse(sorterColumn == 2, 0, 2))
  out = list(firstColumn, secondColumn)
  return(out)
}


# Load sweep data

load("~/Dropbox/oagr_manuscript (1)/WorkDir/SweepFinder/QC_chris/factor_analysis/ray_cutoff_data/SWEEP.LIST.rdata")

mfull = out.list[["Sweeps"]][["0.000625"]]
chr = as.numeric(sapply(strsplit(colnames(mfull), split=":"), function(x) x[1]))
start = sapply(strsplit(colnames(mfull), split=":"), function(x) as.numeric(strsplit(x[2], split="-")[[1]][1])*1e6)
end = sapply(strsplit(colnames(mfull), split=":"), function(x) as.numeric(strsplit(x[2], split="-")[[1]][2])*1e6)

sweepRegions = data.table(chr, start, end)

# Load orignial bed data

setwd("~/Dropbox/oagr_manuscript (1)/WorkDir/new_data/ancient/")

bedData <- BEDMatrix("./ancientStudy.bed")
m = t(as.matrix(bedData))

info = read.table("ancientStudy_arch.ind", stringsAsFactors = F)
group = info$V3
ind = info$V1

snps = fread("ancientStudy.bim")
pos = snps$V4
chromo = snps$V1

ids = gsub(pattern = ".*_(.*$)", replacement = "\\1", colnames(m))

all(ids == ind) # Has to be true if same ordering

# Where are the different notations?
ids[which(ids != ind)]
ind[which(ids != ind)]

# Load Ust Ishim

Ust_bedData <- BEDMatrix("../ustishim/old_HG_studyLoci.bed")
Ust_m = t(as.matrix(Ust_bedData))

Ust_snps = fread("../ustishim/old_HG_studyLoci.bim")
Ust_pos = Ust_snps$V4
Ust_chromo = Ust_snps$V1

any(Ust_pos != pos)  # Both have to be false
any(Ust_chromo != chromo)

# Make sort of pseudo haploid data from diploid calls

apply(Ust_m, 2, table)  # Only WHG_Loschbour and Ust_Ishim_Ust_Ishim have genotype data, others are pseudohaploid!!

Ust_m_allph = Ust_m[,c("MA1_MA1", "CHG_KK1", "CHG_SATP", "Mota_mota", "Clovis_Anzick")]
Ust_m_alldiploid = Ust_m[,c("WHG_Loschbour", "Ust_Ishim_Ust_Ishim")]

Ust_m_alldiploid_pseudohaplodized = apply(Ust_m_alldiploid, 2, function(x) diploid2Pseudohaploid(x))
Ust_m_alldiploid_pseudohaplodized = do.call(cbind, lapply(Ust_m_alldiploid_pseudohaplodized, function(x) do.call(cbind,x)))
colnames(Ust_m_alldiploid_pseudohaplodized) = paste0(rep(colnames(Ust_m_alldiploid), each=2), "_", 1:2)

Ust_m__together_all_pseudohaploid = cbind(Ust_m_allph, Ust_m_alldiploid_pseudohaplodized)


# Correct plarization

wrongPolID = snps$V5 != Ust_snps$V5
Ust_m__together_all_pseudohaploid_polarizedCorrectly = Ust_m__together_all_pseudohaploid
Ust_m__together_all_pseudohaploid_polarizedCorrectly[wrongPolID,] = 2-Ust_m__together_all_pseudohaploid_polarizedCorrectly[wrongPolID,]


# Add Ust Ishim and other ancient HG data to other data
# Important: polarize correctly!!

m = cbind(Ust_m__together_all_pseudohaploid_polarizedCorrectly, m)
group = c(colnames(Ust_m__together_all_pseudohaploid_polarizedCorrectly), group)
group = gsub(pattern = "_[1-2]$", "", group)
ind = c(colnames(Ust_m__together_all_pseudohaploid_polarizedCorrectly), ind)
ids = c(colnames(Ust_m__together_all_pseudohaploid_polarizedCorrectly), ids)


setwd("~/Dropbox/oagr_manuscript (1)/WorkDir/SweepFinder/QC_chris/div_over_time/")

unique(group)



# Load SweepFinder data

fn = dir(path = "../../SweepFinder_OUT", pattern = "ancientStudy_arch.*out", full.names = T, recursive = T)

# Remove populations with little data, i.e. Levant and CentralEurope_HG

#ID1 = !grepl(".*CentralEurope_HG.*", fn)
#ID2 = !grepl(".*Levant.*", fn)

#fn = fn[ID1 & ID2]
fn_plain = sapply(fn, function(x) tail(strsplit(x, split="/")[[1]], n = 1), USE.NAMES = FALSE)
groupSF = gsub(fn_plain, pattern = ".+[sS]tudy_.+?_(.*)_[0-9]+\\.[sS][fF]2out", replacement = "\\1")
groupingSF = gsub(fn_plain, pattern = ".+[sS]tudy_(.+?)_.*_[0-9]+\\.[sS][fF]2out", replacement = "\\1")
chromoSF = gsub(fn_plain, pattern = ".+[sS]tudy_.+?_.*_([0-9]+)\\.[sS][fF]2out", replacement = "\\1")

sfout_files = mclapply((1:length(fn)), function(i) {
  try({cbind(fread(input = fn[i], header = T, sep = "\t"), chr = chromoSF[i], group = groupSF[i], grouping = groupingSF[i])})
}, mc.cores = 8)

sfout_files = do.call(rbind, sfout_files)

sfout_files[, pop := group]







# Look at pseudo haplotype plot -- sweeps at all 40 sweep regions

sortedGroups = c("MA1_MA1", "CHG_KK1", "CHG_SATP", "Mota_mota", "Clovis_Anzick", "WHG_Loschbour", "Ust_Ishim_Ust_Ishim", "SHG", "EHG", "Balkan_HG", "WHG", "EasternEurope_Meso", "Anatolia_Neolithic", "Balkan_Neo", "CentralEurope_Neo", "EasternEurope_Neo", "Iran_Chal_Neo", "WestEuropean_Neo", "UK_Neo", "Scandinavian_Farmers", "CentralEurope_CA", "WestEuropean_CA", "Steppe", "Armenia_BronzeA", "Balkan_BronzeA", "CentralEurope_BronzeA", "EasternEurope_BronzeA", "UK_BronzeA", "WestEuropean_BronzeA")

for (i in 1:nrow(sweepRegions)) {
  
  ichrom = sweepRegions$chr[i]
  istart = sweepRegions$start[i]
  iend = sweepRegions$end[i]
  paddingFactor = 3 # add minus and plus three times region length at both ends
  regionLength = iend - istart
  
  ID = chromo==ichrom & pos > istart - paddingFactor*regionLength & pos < iend + paddingFactor*regionLength
  # groupID = !group %in% c("Levant", "CentralEurope_HG") & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
  groupID = !group %in% c("Levant", "CentralEurope_HG") & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
  
  imageFN = paste0("./all_sweeps_haploImages/haploImage_sweep", i, "_chr", ichrom, "_start", istart/1e6, "_end", iend/1e6, ".jpg")
  jpeg(imageFN, width = 10000, height = 35000, res = 1600)
  par(mar=c(2,0,2,0))
  layout(as.matrix(1:length(unique(sortedGroups))), heights=sapply(sortedGroups, function(x) sum(group==x & groupID)+50))
  for (j in sortedGroups) {
    if (sum(group==j & groupID) == 1) {
      image(z = cbind(m[ID, group==j & groupID], m[ID, group==j & groupID]), x=pos[ID], main=j, col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
      next
    } 
    image(z = m[ID, group==j & groupID], x=pos[ID], main=j, col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
    #image(z = m[ID, group==j & groupID], main=j, col=c("yellow", "black", "black", "black"), cex.main=15)
  }
  dev.off()
  
}



sortedGroups = c("SHG", "EHG", "Balkan_HG", "WHG", "EasternEurope_Meso", "Anatolia_Neolithic", "Balkan_Neo", "CentralEurope_Neo", "WestEuropean_Neo", "UK_Neo", "Steppe", "CentralEurope_BronzeA", "UK_BronzeA")

for (i in 1:nrow(sweepRegions)) {
  
  ichrom = sweepRegions$chr[i]
  istart = sweepRegions$start[i]
  iend = sweepRegions$end[i]
  paddingFactor = 2 # add minus and plus three times region length at both ends
  regionLength = iend - istart
  
  ID = chromo==ichrom & pos > istart - paddingFactor*regionLength & pos < iend + paddingFactor*regionLength
  # groupID = !group %in% c("Levant", "CentralEurope_HG") & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
  groupID = !group %in% c("Levant", "CentralEurope_HG") & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
  
  imageFN = paste0("./all_sweeps_haploImages_JUST_SOME_POPS/haploImage_sweep", i, "_chr", ichrom, "_start", istart/1e6, "_end", iend/1e6, ".jpg")
  jpeg(imageFN, width = 10000, height = 30000, res = 1600)
  par(mar=c(2,0,2,0))
  layout(as.matrix(1:length(unique(sortedGroups))), heights=sapply(sortedGroups, function(x) sum(group==x & groupID)+50))
  for (j in sortedGroups) {
    if (sum(group==j & groupID) == 1) {
      image(z = cbind(m[ID, group==j & groupID], m[ID, group==j & groupID]), x=pos[ID], main=j, col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
      next
    } 
    image(z = m[ID, group==j & groupID], x=pos[ID], main=j, col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
    #image(z = m[ID, group==j & groupID], main=j, col=c("yellow", "black", "black", "black"), cex.main=15)
  }
  dev.off()
  
}


# Plot haplotypes with SF data for some groups/sweeps


i=33
ichrom = sweepRegions$chr[i]
istart = sweepRegions$start[i]
iend = sweepRegions$end[i]
paddingFactor = 1 # add minus and plus three times region length at both ends
regionLength = iend - istart

ID = chromo==ichrom & pos > istart - paddingFactor*regionLength & pos < iend + paddingFactor*regionLength

populID = "Anatolia_Neolithic"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr6_31.6e6_32.6e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 250))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()

populID = "CentralEurope_Neo"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr6_31.6e6_32.6e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 250))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()

populID = "UK_Neo"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr6_31.6e6_32.6e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 250))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()

populID = "UK_BronzeA"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr6_31.6e6_32.6e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 250))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()





i=1
ichrom = sweepRegions$chr[i]
istart = 204.7e6
iend = 204.9e6
paddingFactor = 1 # add minus and plus three times region length at both ends
regionLength = iend - istart

ID = chromo==ichrom & pos > istart - paddingFactor*regionLength & pos < iend + paddingFactor*regionLength

populID = "CentralEurope_Neo"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr1_204.8e6_205.0e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 100))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()

populID = "UK_Neo"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr1_204.8e6_205.0e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 100))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()

populID = "UK_BronzeA"
IDSF = sfout_files$pop == populID & sfout_files$chr == ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
groupID = group %in% c(populID) & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
imageFN = paste0("./EXAMPLES/chr1_204.8e6_205.0e6_", populID, "__SF_HaploPlot.jpeg")
jpeg(imageFN, width = 30000/3, height = 10000/3, res = 1000)
layout(as.matrix(1:2), heights=c(1,1))
par(mar=c(0,4.1,2.1,2.1))
plot(sfout_files[IDSF]$location, sfout_files[IDSF]$LR, type="l", lwd=2, main=populID, xlim=summary(sfout_files[IDSF]$location)[c(1,6)], ylab="", xaxs = "i", xaxt='n', ylim=c(0, 100))
par(mar=c(2,4.1,0,2.1))
image(z = m[ID, groupID], x=pos[ID], col=c("yellow", "black", "black", "black"), cex.main=1, cex.axis=1)
dev.off()



# Diversity function

compute_running_average = function(chr, pos, statistic, windowsize, stepsize) {
  if (windowsize %% stepsize != 0) stop("Windowsize must be multiple of stepsize!")
  data = data.table(chr, pos, statistic)
  setkey(data, chr)
  dataOut = data.table()
  for (chromo in data[, unique(chr)]) {
    seqStarts = seq(0, windowsize-stepsize, stepsize)
    dataOut = rbind(dataOut, do.call(rbind, lapply(seqStarts, function(x) {
      dataChr = data[chr == chromo]
      setkey(dataChr, pos)
      returnValue = dataChr[,pos_cut := cut(pos, breaks = seq(x, max(pos), windowsize), include.lowest = F, dig.lab=10), nomatch=0L][, .(ave = mean(statistic, na.rm=T), min=min(statistic, na.rm=T), max=max(statistic, na.rm=T), median=median(statistic, na.rm=T)), pos_cut, nomatch=0L]
      returnValue[, start := as.numeric(sub(pattern = ".(\\d+),(\\d+).", replacement = "\\1", x = as.character(pos_cut)))+1]
      returnValue[, end := as.numeric(sub(pattern = ".(\\d+),(\\d+).", replacement = "\\2", x = as.character(pos_cut)))]
      returnValue[, midpoint := (start+end)/2]
      returnValue[, chr := chromo]
      return(returnValue)
    }))[order(start)])
  }
  return(dataOut[!is.na(pos_cut)])
}

compNucDiv = function(matrix, pos, chromo, windowsize, stepsize) { 
  freq = rowSums(matrix, na.rm=T)
  sampleSize = rowSums(!is.na(matrix))*2
  perSNPnuclDiv = 2*freq*(sampleSize-freq)/(sampleSize*(sampleSize-1))
  res = compute_running_average(chromo, pos, perSNPnuclDiv, windowsize, stepsize)
  resSampleSize = compute_running_average(chromo, pos, sampleSize, windowsize, stepsize)
  return(cbind(res, ave_sampleSize = resSampleSize$ave))
}

compPolym = function(matrix, pos, chromo, windowsize, stepsize) { 
  freq = rowSums(matrix, na.rm=T)
  sampleSize = rowSums(!is.na(matrix))*2
  perSNPpolym = 1-as.numeric(freq==sampleSize | freq==0)
  res = compute_running_average(chromo, pos, perSNPpolym, windowsize, stepsize)
  resSampleSize = compute_running_average(chromo, pos, sampleSize, windowsize, stepsize)
  return(cbind(res, ave_sampleSize = resSampleSize$ave))
}



# Look at diversity at this sweep

for (i in 1:nrow(sweepRegions)) {
  print(i)
  ichrom = sweepRegions$chr[i]
  istart = sweepRegions$start[i]
  iend = sweepRegions$end[i]
  #paddingFactor = 3 # add minus and plus three times region length at both ends
  #regionLength = iend - istart
  paddingFactor = 1
  regionLength = 2e6
  
  ID = chromo==ichrom & pos > istart - paddingFactor*regionLength & pos < iend + paddingFactor*regionLength
  # groupID = !group %in% c("Levant", "CentralEurope_HG") & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 0.5) # Less than 50% missing data
  groupID = !group %in% c("Levant", "CentralEurope_HG") & (colSums(is.na(m[ID, ]))/nrow(m[ID, ]) < 1) # Less than 100% missing data --all data
  
  IDsf = as.numeric(sfout_files$chr)==ichrom & sfout_files$location > istart - paddingFactor*regionLength & sfout_files$location < iend + paddingFactor*regionLength
  
  nuclDiv = data.table()
  polym = data.table()
  for (j in names(table(group[groupID]))) {
    nuclDiv = rbind(nuclDiv, cbind(pop=j, compNucDiv(m[ID, group==j & groupID], pos[ID], chromo[ID], 0.5e6, 0.01e6)))
    polym = rbind(polym, cbind(pop=j, compPolym(m[ID, group==j & groupID], pos[ID], chromo[ID], 0.5e6, 0.01e6)))
  }

  nuclDiv$pop_F = factor(nuclDiv$pop, levels = c("SHG", "SEHG", "EHG", "Balkan_HG", "WHG", "Steppe", "EasternEurope_Meso", "Anatolia_Neolithic", "Balkan_Neo", "CentralEurope_Neo", "EasternEurope_Neo", "Iran_Chal_Neo", "WestEuropean_Neo", "UK_Neo", "Scandinavian_Farmers", "CentralEurope_CA", "WestEuropean_CA", "Armenia_BronzeA", "Balkan_BronzeA", "CentralEurope_BronzeA", "EasternEurope_BronzeA", "UK_BronzeA", "WestEuropean_BronzeA"))
  sfout_files$pop_F = factor(sfout_files$pop, levels = c("SHG", "SEHG", "EHG", "Balkan_HG", "WHG", "Steppe", "EasternEurope_Meso", "Anatolia_Neolithic", "Balkan_Neo", "CentralEurope_Neo", "EasternEurope_Neo", "Iran_Chal_Neo", "WestEuropean_Neo", "UK_Neo", "Scandinavian_Farmers", "CentralEurope_CA", "WestEuropean_CA", "Armenia_BronzeA", "Balkan_BronzeA", "CentralEurope_BronzeA", "EasternEurope_BronzeA", "UK_BronzeA", "WestEuropean_BronzeA"))

  nuclDiv[,max_ave_sampleSize := max(ave_sampleSize), pop]
  polym[,max_ave_sampleSize := max(ave_sampleSize), pop]

  highSampPops = c("EHG", "Balkan_HG", "Steppe", "Anatolia_Neolithic", "Balkan_Neo", "CentralEurope_Neo", "UK_Neo", "CentralEurope_CA", "WestEuropean_CA", "UK_BronzeA", "CentralEurope_BronzeA", "Balkan_BronzeA")

  #ggplot(data=polym) + geom_line(aes(midpoint, ave)) + theme_bw() + facet_wrap("pop") + geom_vline(aes(xintercept=32e6), lty=2, col="grey") + 
  #  geom_line(aes(midpoint, ave_sampleSize/max_ave_sampleSize), col="black")
  
  imageFN_div = paste0("./all_sweeps_diversity/diversity_sweep", i, "_chr", ichrom, "_start", istart/1e6, "_end", iend/1e6, ".jpg")
  
  jpeg(file=imageFN_div, width = 300, height = 1500)
    maxLRvalue = sfout_files[IDsf & pop %in% highSampPops][LR==max(LR), location]
    pl = ggplot(data=nuclDiv[pop %in% highSampPops]) + geom_line(aes(midpoint, ave)) + theme_bw() + facet_grid(pop_F~.) + geom_vline(aes(xintercept=maxLRvalue), lty=2, col="darkgrey") + 
      geom_line(data = sfout_files[IDsf & pop %in% highSampPops], aes(location, LR/max(LR)*0.5), col="red") + xlab("Position (chr 6)")
    print(pl)
  dev.off()

  #ggplot(data=sfout_files[IDsf]) + geom_line(aes(location, LR, col=pop)) + theme_bw() + facet_wrap("pop") + geom_vline(aes(xintercept=32e6), lty=2, col="grey") 

}

