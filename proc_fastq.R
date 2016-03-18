#########################
## FASTQ Quality Plots ##
#########################
## Author: Li Pidong
## Date: 2016-03-17


# process data ------
#' Processing of single fastq file
#'
#' @param fastq_file
#'
#' @return list
#' quality_score: the matrix of  quality score
#' bases : the matrix of bases
#' nreads : the number of reads in this file
proc_fastq <- function(fastq_file){
  library(ShortRead)
  library(Biostrings)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  fq <- readFastq(fastq_file)
  nReads <- length(fq) # Total number of reads in fastq file
  ## If reads are not of constant width then inject them into a matrix pre-populated with
  if(length(unique(width(fq))) == 1) {
    q <- as.matrix(PhredQuality(quality(fq)))# qs矩阵
    s <- as.matrix(sread(fq))# 碱基矩阵
  } else {# 将不足的碱基和质量得分补全
    mymin <- min(width(fq));
    mymax <- max(width(fq))
    s <- matrix("N", length(fq), mymax)
    q <- matrix(NA, length(fq), mymax)
    for(i in mymin:mymax) {
      index <- width(fq)==i
      if(any(index)) {
        s[index, 1:i] <- as.matrix(DNAStringSet(sread(fq)[index], start=1, end=i))
        q[index, 1:i] <- as.matrix(PhredQuality(quality(fq))[index])
      }
    }
  }
  return(list(quality_score=q, bases=s, nreads=nReads, seq_length=ncol(q)))
}


# plot ------

#' point plot distribution of quality score
#'
#' @param fastq_data  from the function of proc_fastq
#'
#' @return point plot
point_plot <- function(fastq_data){
  qs <- fastq_data$quality_score
  qs_reshape <- melt(qs)
  colnames(qs_reshape) <- c('reads', 'pos', 'quality_score')
  p <- ggplot(data=qs_reshape, aes(pos, quality_score) )
  plot_out <- p+geom_point(alpha=0.1, size=3, fill=rgb(15, 115, 175,maxColorValue = 255),color=rgb(15, 115, 175,maxColorValue = 255),
                     shape=19)+
    xlab('Position along reads')+
    ylab('Quality')+
    ggtitle('Distribution of qualities')+
    theme_bw()+
    theme(plot.title=element_text(vjust=1))

}

#' Title calulate the base freqency
#'
#' @param fastq_data from the function of proc_fastq
#'
#' @return base frequency
base_freq <- function(fastq_data){
  bases <- fastq_data$bases
  bstats <- apply(bases, 2, function(x) table(factor(x, levels=c("A", "C", "G", "T","N"))))
  colnames(bstats) <- 1:length(bstats[1,])
  # the percentage of per base
  bstats <- t(apply(bstats, 1, function(x) x/colSums(bstats)))
  bstats <- melt(bstats)
  colnames(bstats) <- c('base', 'pos', 'freq')
  bstats$pos <- factor(bstats$pos, levels=unique(bstats$pos), ordered=T)
  bstats
}

#'Title per base frequency plot------
#'
#' @param fastq_data from the function of proc_fastq
#'
#' @return line plot
base_freq_plot <- function(fastq_data){

  bstats <- base_freq(fastq_data)
  out <- ggplot(bstats, aes(pos, freq, group=base, color=base)) +
    scale_x_discrete(breaks=c(1, seq(0, length(unique(bstats$pos)), by=10)[-1])) +
    geom_line(stat="identity") +
    ylim(0,0.5)+
    xlab('Position along reads')+
    ylab("Proportion")+
    scale_color_discrete(name="Bases")+
    theme_bw()+
    ggtitle('Base percentage composition along reads')+
    theme(plot.title=element_text(vjust=1))
}

#'Title box_plot of the distribution of quality score
#'
#' @param fastq_data from the function of proc_fastq
#'
#' @return ggplot
qs_box_plot <- function(fastq_data){
  quality_score <- fastq_data$quality_score
  row.names(quality_score) <- paste("s", 1:length(quality_score[,1]), sep="");
  colnames(quality_score) <- 1:length(quality_score[1,])
  bpl <- boxplot(quality_score, plot=FALSE)# 获取数值
  astats <- data.frame(bpl$names, t(matrix(bpl$stats, dim(bpl$stats))))
  colnames(astats) <- c("pos", "min", "low", "mid", "top", "max")
  astats[,1] <- factor(astats[,1], levels=unique(astats[,1]), ordered=TRUE)
  # plot
  out <- ggplot(astats, aes(x=pos, ymin = min, lower = low, middle = mid, upper = top, ymax = max)) +
    geom_boxplot(stat = "identity",
                 color="#999999",
                 fill="#56B4E9",
                 varwidth = T) +
    scale_x_discrete(breaks=c(1, seq(0, length(astats[,1]), by=10)[-1])) +
    theme(legend.position = "none", plot.title = element_text(size = 12))+
    theme_bw()+
    xlab('Position along reads')+
    ylab('Quality')+
    ggtitle('Distribution of qualities')+
    theme(plot.title=element_text(vjust=1))
}

#' Title per position  base freqency bar plot
#'
#' @param fastq_data  from the function of proc_fastq
#'
#' @return bar plot
base_freq_bar_plot <- function(fastq_data){
    bstats <- base_freq(fastq_data)
    out <- ggplot(bstats, aes(pos, freq, fill=base), color="black") +
    scale_x_discrete(breaks=c(1, seq(0, length(unique(bstats$pos)), by=10)[-1])) +
    geom_bar(stat="identity") +
    # scale_y_continuous(expand = c(0,0))+
    ylab("Proportion")+
    xlab('Position along reads')+
    ggtitle('Base percentage composition along reads')+
    theme_bw()+
    theme(plot.title=element_text(vjust=1))
}


#' Title Per position average quality score of each base type
#'
#' @param fastq_data from the function of proc_fastq
#'
#' @return bar plot
base_qs_per_pos_plot <- function(fastq_data){
  q <- fastq_data$quality_score
  s <- fastq_data$bases
  A <- q; A[s %in% c("T", "G", "C", "N")] <- NA; A <- colMeans(A, na.rm=TRUE)
  T <- q; T[s %in% c("A", "G", "C", "N")] <- NA; T <- colMeans(T, na.rm=TRUE)
  G <- q; G[s %in% c("T", "A", "C", "N")] <- NA; G <- colMeans(G, na.rm=TRUE)
  C <- q; C[s %in% c("T", "G", "A", "N")] <- NA; C <- colMeans(C, na.rm=TRUE)
  cstats <- data.frame(qs=c(A, C, G, T), bases=rep(c("A", "C", "G", "T"), each=length(A)), pos=rep(1:length(A),4))

  out <- ggplot(cstats, aes(pos, qs, group=bases, color=bases)) +
    geom_line() +
    scale_x_discrete(breaks=c(1, seq(0, length(unique(cstats$pos)), by=10)[-1]))+
    xlab('Position along reads')+
    ylab('Quality score')+
    ggtitle('Per position average quality of each base type')+
    theme_bw()+
    theme(plot.title=element_text(vjust=1))
}

#' Title average quality score per read
#'
#' @param fastq_data
#'
#' @return ggplot
per_read_avg_qs_plot <- function(fastq_data){
  # average quality score per read----
  qs <- fastq_data$quality_score
  qs_reshape <- melt(qs)
  colnames(qs_reshape) <- c('reads', 'pos', 'quality_score')
  qs_avg <- qs_reshape %>%
    group_by(reads) %>%
    summarise(avg_per_read = mean(quality_score))
  out <- ggplot(qs_avg, aes(avg_per_read)) +
    # geom_density(stat = 'bin', color='#DC0000')+
    geom_histogram(aes(y=..density..), fill="#0072B2", stat="bin", color='white', binwidth=1.5)+
    theme_bw()+
    ylab("Density")+
    xlab('The mean of  read quality')+
    ggtitle('Quality score distribution over all reads')+
    theme_bw()+
    theme(plot.title=element_text(vjust=1))+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0), limits=c(0,0.16))
}

#' GC content distribution
#'
#' @param fastq_data from the function of proc_fastq
#'
#' @return ggplot
gc_plot <- function(fastq_data){
  bases <- fastq_data$bases
  nbase <- ncol(bases)
  dstats <- sapply(1:nrow(bases), function(ii){
    sum(bases[ii,] %in% c('C','G', 'c', 'g'))/nbase
  })
  dstats <- data.frame(gc_content=dstats)
  nrom_data <- data.frame(x=rnorm(nrow(dstats),mean=mean(dstats$gc_content), sd=0.08))
  out <- ggplot(dstats) +
    geom_density(aes(x=gc_content, y=..density..), stat = 'bin')+
    stat_function(fun=dnorm,aes(x), args=list(mean=mean(dstats$gc_content), sd=0.08), data=nrom_data, color='red')+
    theme_bw()+
    theme(axis.text.y = element_blank())+
    xlab('Mean GC content')+
    ggtitle('GC content distribution over all reads')+
    theme(plot.title=element_text(vjust=1))
}

#' point plot distribution of quality score
#'
#' @param fastq_data
#'
#' @return point plot
point_plot_1 <- function(fastq_data){
  qs <- fastq_data$quality_score
  qs_reshape <- melt(qs)
  colnames(qs_reshape) <- c('reads', 'pos', 'quality_score')
  plot(qs_reshape$pos, qs_reshape$quality_score)
  # p <- ggplot(data=qs_reshape, aes(pos, quality_score) )
  # plot_out <- p+geom_point(alpha=0.1, size=3, fill=rgb(15, 115, 175,maxColorValue = 255),color=rgb(15, 115, 175,maxColorValue = 255),
  #                          shape=19)+
  #   xlab('Position along reads')+
  #   ylab('Quality')+
  #   ggtitle('Distribution of qualities')+
  #   theme_bw()+
  #   theme(plot.title=element_text(vjust=1))
  #
}

