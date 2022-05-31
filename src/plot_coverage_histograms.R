#!/usr/bin/env Rscript

library(optparse)

# read user info
option_list = list(
  make_option(c("-c", "--cvgFile"), type="character", default=NA, help="coverage file from bedtools"),
  make_option(c("-b", "--bedFile"),  type="character", default=NA, help="input bed file", metavar="character"),
  make_option(c("-o", "--outputBasename"),type="character", default=NA, help="basename for plots created", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# test if there is at least one argument: if not, return an error
if (length(option_list)!=3) {
  print_help(opt_parser)
  stop("", call.=FALSE)
}

library(ggpubr)
library(ggplot2)
library(tidytext)

#test
options(bitmapType='cairo')

#read bed file
bed<-read.delim(file=opt$bedFile, sep="\t",header=F)

colnames(bed)<-c("chrom","start","stop","pool")
bed<-bed[order(bed$pool,bed$chrom,bed$start),]
rownames(bed)<-paste(bed$chrom,":",bed$start,"-",bed$stop,sep="")

df.all<-NULL
id<-opt$outputBasename


#read coverage histogram
hist<-read.delim(file=opt$cvgFile,sep="\t",as.is=T,header=F)

colnames(hist)<-c("chrom","start","stop","pool","depth","bases","size","proportion")
hist<-hist[hist$chrom != "all",]
hist$interval<-paste(hist$chrom,":",hist$start,"-",hist$stop,sep="")

#calculate mean coverage metric
mean_coverage<-by(hist,hist$interval,function(x){sum(x$depth*x$bases)/as.numeric(x$size[1])})
intervals<-as.vector(names(mean_coverage))
mean_coverage<-as.vector(mean_coverage)

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

### set intreval as factor to order by intervsla in the bed file
df.all$interval<-factor(df.all$interval,levels=rownames(bed))
df.all$pool<-bed[df.all$interval,]$pool

#save.image(file="image.RData")
############# All plot
g0<-ggplot(df.all[df.all$metric=="cvg_mean",], aes(x=interval,y=value,col=pool)) +
  geom_bar(stat="identity") +
  theme( axis.text.x = element_blank()) +
  labs(title=paste("                                                                                                  ",id, "\nMean Interval coverage", sep = "")) + xlab("interval") + ylab("depth")
ggsave(g0,file=paste0(id,"_mean_interval_coverage.png"),dev="png",height=10,width=15)

############# Percent coverage
cvg0<-df.all[df.all$metric=="cvg0",]
cvg_mean<-df.all[df.all$metric=="cvg_mean",]
percent_intervals_with_coverage<-aggregate(value~id + pool, data=cvg_mean,function(x){length(x[x>0])/length(x)*100})

g1<-ggplot(percent_intervals_with_coverage,aes(y=value,x=pool,col=pool)) +
  geom_bar(stat="identity") +
  #facet_wrap(~id,ncol=3) +
  theme(axis.text.x = element_blank()) +
  labs(title=paste ("Percent of Intervals with coverage", sep = "")) + xlab("subset") + ylab("percent")
ggsave(g1,file=paste0(id,"_percent_intervals_w_coverage.png"),dev="png",height=10,width=15)

############# Sorted coverage
g2<-ggplot(df.all[df.all$metric == "cvg_mean",],aes(x=reorder_within(interval,value,list(id)),y=value)) + geom_point() +
  scale_y_log10()+
  theme(axis.text.x = element_blank()) +
  guides(x = "none") +
  labs(title=paste ("Mean interval coverage sorted", sep = "")) + xlab("interval") + ylab("depth")

for (pool in unique(df.all$pool)){
  if (length(which(df.all$pool == pool)) > 10000 ){
    #print(pool)
    df2<-df.all
    df2$set<-ifelse(df2$pool== pool, pool,"Other pools")

    subset1<-rownames(bed[bed$pool!=pool,])
    subset2<-sample(rownames(bed[bed$pool==pool,]),4000)
    subset<-c(subset1,subset2)

    df2<-df2[df2$interval %in% subset,]

    g2<-ggplot(df2[df2$metric == "cvg_mean",],aes(x=reorder_within(interval,value,list(id,set)),y=value)) + geom_point() +
      facet_wrap(~set,ncol=2,scales="free_x") +
      scale_y_log10()+
      theme(axis.text.x = element_blank()) +
      guides(x = "none") +
      labs(title=paste ("Mean interval coverage sorted", sep = "")) + xlab("interval") + ylab("depth")

    ############# PLOT LARGE POOL
    large_set <- rownames(bed[bed$pool== pool,])
    df.largePool <- df.all[df.all$interval %in% large_set,]

    g3<-ggplot(df.largePool[df.largePool$metric=="cvg_mean",], aes(x=interval,y=value,col=pool)) +
      geom_bar(stat="identity") +
      theme(axis.text.x = element_blank()) +
      labs(title=paste("                                                                                                  ",id, "\nMean interval coverage - Large pool", sep = "")) + xlab("interval") + ylab("depth")
    ggsave(g3,file=paste0(id,"_mean_interval_coverage_large_pool.png"),dev="png",height=10,width=15)


    ############# Pool and Subsampled Large pool
    set1<-rownames(bed[bed$pool!=pool,])
    set2<-sample(rownames(bed[bed$pool==pool,]),800)
    set<-c(set1,set2)
    df.set<-df.all[df.all$interval %in% set,]

    ### SAMPLE FURTHER
    g4<-ggplot(df.set[df.set$metric=="cvg_mean",],aes(x=interval,y=value,col=pool)) +
      geom_bar(stat="identity") +
      #facet_wrap(~id,ncol=4) +
      theme(axis.text.x = element_blank()) +
      labs(title=paste ("Mean interval coverage - Pools and subsampled large pool", sep = "")) + xlab("interval") + ylab("depth")
    ggsave(g4,file=paste0(id,"_mean_interval_coverage_pools_subsampled_largepool.png"),dev="png",height=10,width=15)
  }
}

ggsave(g2,file=paste0(id, "_mean_interval_coverage_sorted.png"),dev="png",height=10,width=15)
