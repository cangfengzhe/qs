#########################
## FASTQ Quality Plots ##
#########################
## Author: Li Pidong
## Date: 2016-03-17

library(ShortRead)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(reshape2)
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
  return(list(quality_score=q, bases=s, nreads=nReads))
}


# plot ------

#' point plot distribution of quality score
#'
#' @param quality_score the matrix of quality_score from the function of proc_fastq
#'
#' @return point plot
point_plot <- function(quality_score){
  qs <- quality_score
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

#'title per base percentage plot------
#'
#' @param bases the matrix of bases from the function of proc_fastq
#'
#' @return ggplot
base_percentage_plot <- function(bases){
  bstats <- apply(bases, 2, function(x) table(factor(x, levels=c("A", "C", "G", "T","N"))))
  colnames(bstats) <- 1:length(bstats[1,])
# the percentage of per base
  bstats <- t(apply(bstats, 1, function(x) x/colSums(bstats)))
  bstats <- melt(bstats)
  colnames(bstats) <- c('base', 'pos', 'freq')
  bstats$pos <- factor(bstats$pos, levels=unique(bstats$pos), ordered=T)

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

#' box_plot of the distribution of quality score
#'
#' @param quality_score the matrix of quality_score from the function of proc_fastq
#'
#' @return ggplot
box_plot <- function(quality_score){
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
    scale_x_discrete(breaks=c(1, seq(0, length(fq[['qs_box']][,1]), by=10)[-1])) +
    theme(legend.position = "none", plot.title = element_text(size = 12))+
    theme_bw()+
    xlab('Position along reads')+
    ylab('Quality')+
    ggtitle('Distribution of qualities')+
    theme(plot.title=element_text(vjust=1))
}

#' Title the bar plot of  base percentage
#'
#' @param quality_score the matrix of quality_score from the function of proc_fastq
#' @param bases the matrix of bases from the function of proc_fastq
#'
#' @return ggplot
base_percentage_bar_plot <- function(quality_score, bases){
  q <- quality_score
  s <- bases
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

#' GC content distribution
#'
#' @param bases the matrix of bases from the function of proc_fastq
#'
#' @return
#' @export
#'
#' @examples
gc_plot <- function(bases){
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
#' @param quality_score the matrix of quality_score from the function of proc_fastq
#'
#' @return point plot
point_plot_1 <- function(quality_score){
  qs <- quality_score
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

