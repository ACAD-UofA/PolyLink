library(stringi)
library(data.table)
library(matrixStats)
library(foreach)
library(ggplot2)
library(RColorBrewer)

#set working directory
work.dir <- "~/Dropbox/Selection_test/"
setwd(work.dir)

#path to raw files (genes + pathway info)
raw.path <- file.path(work.dir, "Raw_files")

#PolyLink output directory
out.path <- file.path(work.dir, "PolyLink/Output")

ll <- list.files(out.path)

#enter all PolyLink output data
d.newmeth <- foreach(j=ll, .combine=rbind) %do% {
   fp <- file.path(out.path, j, "results")
   ll2 <- list.files(fp)
   f.now <- ll2[grep("setscore", ll2)]
   if(length(f.now)>0){
      ff <- fread(file.path(fp, f.now))
      data.table(Pop=j, ff)
   }
}

#add source info

go.cats <- fread(file.path(raw.path, "PolyLink_GO_pathways_HG19.txt"))
d.newmeth <- merge(d.newmeth, go.cats[!duplicated(Pathway), .(Pathway, Source)], 
                   by.x="setName", by.y="Pathway")

#add effective population size
d.newmeth <- merge(n.eff, d.newmeth, by.x="Population", by.y="Pop")

#add new pop'n label combining eff N
d.newmeth[, Pop.name:=paste0(Population, " [", round(eff.N, 2), "]")]
ord.eff.n <- d.newmeth[!duplicated(Population)][order(eff.N), 
                                                paste0(Population, " [", round(eff.N, 2), "]")]
d.newmeth[, Pop.name:=factor(Pop.name, levels=ord.eff.n)]

#redo q value correction, omitting Enard viral categories
d.newmeth[, setQ.OLD:=setQ]
d.newmeth[!(Source %in% "Enard_Viral"), setQ:=qvalue(setP)$qvalues, by="Population"]

##---------------------------------------------------------------------------##
#  diagnostic plots showing top setScore distributions across populations
##---------------------------------------------------------------------------##

plot.path <- file.path(work.dir, "Plots")

#plot scores vs # genes
d.newmeth[, p.sig:=cut(-log10(setP), c(0,1,2,3,4,5,Inf),
                             c("[0,1]", "(1,2]", "(2,3]",
                              "(3,4]", "(4,5]", ">5"), include.lowest=T)]
d.newmeth[, q.sig:=factor(setQ<0.05)]

#number of columns in plots (alter to preferences)
n <- d.newmeth[, length(unique(Population))]


ggplot(d.newmeth, aes(x=N, y=setScore, col=p.sig)) +
   geom_point(size=0.25) + scale_x_log10() +
   facet_wrap(~Pop.name, ncol=n, scale="free_y") +
   scale_color_manual(breaks=c("[0,1]", "(1,2]", "(2,3]",
                              "(3,4]", "(4,5]", ">5"),
                        values = c("grey50", brewer.pal(5, "YlOrRd")))+
   xlab("# genes in pathway") + ylab("Sumstat") + theme_bw() +
   theme(axis.text.x=element_text(size=6, hjust=1),
         axis.text.y=element_text(size=6),
         strip.text=element_text(size=8),
         panel.border=element_blank(),
         panel.spacing.x=unit(0.025, "lines"),
         panel.spacing.y=unit(0.01, "lines")) +
   guides(color=guide_legend(override.aes=list(size=2), title="P Value"))

ggsave(filename=file.path(plot.path, "setScore_by_setSize_by_setP.genes.pdf"),
    width=8, height=6, dpi=300)

#p value distribution
ggplot(d.newmeth, aes(x=setP)) +
   geom_histogram(aes(y =..density..),
                  breaks=seq(0,1,0.01)) +
   facet_wrap(~Pop.name, ncol=n) +
   xlab("p value") + ylab("Density") +
   theme_bw() +
   theme(axis.text.x=element_text(size=6, hjust=1),
         axis.text.y=element_text(size=6),
         strip.text=element_text(size=8),
         panel.border=element_blank(),
         panel.spacing.x=unit(0.025, "lines"),
         panel.spacing.y=unit(0.01, "lines"))

ggsave(filename=file.path(plot.path, "PolyLink.p.value.distribution.pdf"),
       width=8, height=6, dpi=300)



##---------------------------------------------------------------------------##
#  plots showing top pathways across populations
##---------------------------------------------------------------------------##

small.pops <- d.newmeth[floor(eff.N) <= 10, unique(Population)]
d.newmeth[, Class:="Ancient"]
d.newmeth[Population %in% small.pops, Class:="Small"]

#determine significant pathways (q<0.05)
sig.paths <- d.newmeth[!(Class == "Small" | Population %in% c("CHB", "YRI"))
                       & setQ<0.05]
d <- d.newmeth[setName %in% sig.paths$setName,
               .(Class, Population,
                 setName, p=setP, q=setQ)]

#add q.sig symbol
d[, q.sig.symb:=""]
d[q<=0.1, q.sig.symb:="*"]
d[q<=0.05, q.sig.symb:="**"]
d[q<=0.01, q.sig.symb:="!"]

#add binned p-values
d[, `-log10(p)`:= cut(-log10(p), c(0:5, 100), include.lowest=T,
                      c("[0,1]", "(1,2]", "(2,3]", "(3,4]", "(4,5]", ">5"))]

d[, Q:= cut(q, c(0, 0.01, 0.05, 0.1, 0.2, 0.5, 1), include.lowest=T)]


for(cl in list("", "Small")){

   if(cl=="") d.now <- d[!(Class == "Small")]
   if(cl=="Small") d.now <- d

   ggplot(d.now, aes(x=Population, y=setName))+
      geom_tile(aes(fill=`-log10(p)`), color="black")+
      geom_text(aes(label=q.sig.symb), size=5, col="white")+
      facet_grid(~Class, scales="free", space="free")+
      scale_fill_manual(breaks=c("[0,1]", "(1,2]", "(2,3]",
                                 "(3,4]", "(4,5]", ">5"),
                        #values = c(brewer.pal(6, "Greys")[1:3], brewer.pal(6, "OrRd")[4:6]))+
                        values = c("grey95", brewer.pal(5, "YlOrRd")))+
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      labs(fill="-log10(p)") +
      theme_bw()+
      theme(axis.text.x = element_text(size=7, angle=45, hjust=1),
            axis.text.y=element_text(size=7),
            legend.text=element_text(size=8),
            legend.title=element_text(size=10),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            #panel.border=element_blank(),
            panel.border=element_rect(color="black", fill=NA, size=1),
            strip.background=element_rect(color="black", size=1),
            panel.spacing.x=unit(0.01, "lines"),
            panel.spacing.y=unit(0.01, "lines"),
            strip.text.x=element_text(size=7),
            strip.text.y=element_text(size=7),
            legend.position="bottom",
            legend.key.size=unit(1, "line"),
            legend.margin=margin(t=0, r=0, b=-0.05, l=0, unit="cm")) +
      guides(fill=guide_legend(nrow=1), size=4)

   out <- paste0(plot.path, "/heatmap_PolyLink_SweepFinder", cl, ".pdf")

   ggsave(filename=out, width=8, height=4, dpi=300)
}
