#!/usr/bin/env Rscript

library(optparse)

# read user info
option_list = list(
  make_option(c("-c", "--cvgFile"), type="character", default=NA, help="coverage file from bedtools"),
  make_option(c("-b", "--bedFile"),  type="character", default=NA, help="input bed file", metavar="character"),
  make_option(c("-o", "--outputBasename"),type="character", default=NA, help="basename for plots created", metavar="character"),
  make_option(c("-l", "--listSplit"),type="character", default=NA, help="list of views to split", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# error if missing argument: if not, return an error
if (is.na(opt$cvgFile) || is.na(opt$bedFile) || is.na(opt$outputBasename)) {
  print("Missing argument")
  print_help(opt_parser)
  stop("", call.=FALSE)
}

library(ggpubr)
library(ggplot2)
library(tidytext)
options(bitmapType='cairo')

#read bed file
bed<-read.delim(file=opt$bedFile, sep="\t",header=F)
#read coverage histogram
hist<-read.delim(file=opt$cvgFile,sep="\t",as.is=T,header=F)

if (ncol(bed) == 3 && ncol(hist) == 7) {
  #reorganize and add dummy pool for bed and hist when there are no pools in bedfile
  colnames(bed)<-c("chrom","start","stop")
  bed$pool <- "pool_1"
  bed <- bed[, c("chrom","start","stop","pool")]

  colnames(hist)<-c("chrom","start","stop","depth","bases","size","proportion")
  hist$pool <- "pool_1"
  hist <- hist[, c("chrom","start","stop","pool","depth","bases","size","proportion")]

} else {
  colnames(bed)<-c("chrom","start","stop","pool")
  colnames(hist)<-c("chrom","start","stop","pool","depth","bases","size","proportion")
}


bed<-bed[order(bed$pool,bed$chrom,bed$start),]
rownames(bed)<-paste(bed$chrom,":",bed$start,"-",bed$stop,sep="")

df.all<-NULL
id<-opt$outputBasename

hist<-hist[hist$chrom != "all",]
hist$interval<-paste(hist$chrom,":",hist$start,"-",hist$stop,sep="")

#calculate mean coverage metric
mean_coverage<-by(hist,hist$interval,function(x){sum(x$depth*x$bases)/as.numeric(x$size[1])})
intervals<-as.vector(names(mean_coverage))
mean_coverage<-as.vector(mean_coverage)
max_coverage <- max(mean_coverage)

#create dataframe
df<-data.frame(id=id,interval=intervals,metric="cvg_mean",value=mean_coverage)
df.all<-rbind(df.all,df)

#metric for no coverage
intervals<-unique(hist$interval)
no_coverage<-rep(0,length(intervals))
names(no_coverage)<-intervals
hist.nocoverage<-hist[hist$depth==0,]
no_coverage[hist.nocoverage$interval]<-hist.nocoverage$proportion*100
#add metric to coverage
df<-data.frame(id=id,interval=intervals,no_coverage=no_coverage)
df<-data.frame(id=id,interval=intervals,metric="cvg0",value=no_coverage)
df.all<-rbind(df.all,df)

#CALCULATE PERCENT COVERED
hist$covered<-ifelse(hist$depth>0,1,0)
hist$proportion2<-hist$proportion*hist$covered
percent_covered <- by(hist,hist$interval,function(x){sum(x$proportion2)})
intervals<-as.vector(names(percent_covered))
percent_covered<-as.vector(percent_covered)
df<-data.frame(id=id,interval=intervals,metric="pct_cvd",value=percent_covered)
df.all<-rbind(df.all,df)

### set intreval as factor to order by intervsla in the bed file
df.all$interval<-factor(df.all$interval,levels=rownames(bed))
df.all$pool<-bed[df.all$interval,]$pool

to_plot <- list()

if (!is.na(opt$listSplit)) {
#if (!is.na(test)) {
  list_to_split <- strsplit(opt$listSplit, split = ",")
  #list_to_split <- strsplit(test, split = ",")
  print(list_to_split)
  i = 1

  for (pool in list_to_split) {
    print(pool)
    print(nrow(df.all))
    #create a subset with only the pool
    df.subset <- df.all[df.all$pool == pool,]
    print(nrow(df.subset))

    #remove pool from the df.all
    df.all <- df.all[df.all$pool != pool,]

    #append the pool
    to_plot[[i]] <- df.subset
    i = i + 1
  }

  #append the all pool
  to_plot[[length(list_to_split)+1]] <- df.all
} else {
    to_plot[[1]] <- df.all
}


g0_list <- list()
g1_list <- list()
g2_list <- list()

index = 1

for (df_plot in to_plot) {

  ############# All plot
  g0<-ggplot(df_plot[df_plot$metric=="cvg_mean",], aes(x=interval,y=value,col=pool)) +
    geom_bar(stat="identity") +
    theme( axis.text.x = element_blank()) +
    xlab("interval") + ylab("depth")
  g0_list[[index]] <- g0

  ############# Percent coverage
  g1<-ggplot(df_plot[df_plot$metric=="pct_cvd",], aes(x=as.factor(pool), y=value, col=pool)) +
    geom_boxplot(fill="slateblue", alpha=0.2) +
    geom_jitter() +
    xlab("pool") + ylab("proportion")+
    theme(axis.text.x = element_blank()) #+
  g1_list[[index]] <- g1

  ############# Sorted coverage
  g2<-ggplot(df_plot[df_plot$metric == "cvg_mean",],aes(x=reorder_within(interval,value,list(id)),y=value,col=pool)) +
    #scale_x_continuous(limits = c(0, 10000)) +
    ylim(0, max_coverage) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_blank()) +
    guides(x = "none") +
    xlab("interval") + ylab("depth")
  g2_list[[index]] <- g2

  index = index + 1
}

g0all <-ggarrange(plotlist=g0_list, nrow = length(g0_list))
g0all <- annotate_figure(g0all, top = text_grob("Mean interval coverage ordered by position"))
ggsave(g0all,file=paste0(id, "_mean_interval_coverage.png"),dev="png",height=10,width=15)

g1all <-ggarrange(plotlist=g1_list, nrow = length(g1_list) )
g1all <- annotate_figure(g1all, top = text_grob("Proportion of interval covered"))
ggsave(g1all,file=paste0(id, "_interval_proportion_covered.png"),dev="png",height=10,width=15)

g2all <-ggarrange(plotlist=g2_list, ncol = length(g2_list))
g2all <- annotate_figure(g2all, top = text_grob("Mean interval coverage sorted by depth"))
ggsave(g2all,file=paste0(id, "_mean_interval_coverage_sorted.png"),dev="png",height=10,width=15)
