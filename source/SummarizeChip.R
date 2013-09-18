# summarize ChIP-seq data
SummarizeChip<-function(bamchop, verbose=TRUE,
                        shift=c(seq(0, 250, 10), seq(275, 500, 25), seq(550, 1000, 50)),
                        height=c(10, 25, 50, 100, 200, 500, 1000, 5000, 10000, 100000)) {
# bamchop   A 'bamchop' object
  library(IRanges);
  if (verbose) print('Summarizing ChIP-seq data ...');
  chip<-list();
  
  stats<-bamchop$stats;
  reads<-bamchop$mapping.location;
  genome<-bamchop$genome;  
  
  frag.len<-bamchop$extra.param$fragment.length; # average fragment length, used to extend reads
  
  ################################################## strand-strand correlation ################################################## 
  if (verbose) print('Summarizing strand-strand correlation ...');
  x<-sapply(reads, function(r) StrandCorr(r, resize=1, shift=shift, nodup=TRUE));
  w<-unlist(x[1,]);
  corr<-data.frame(do.call('rbind', x[2,]));
  colnames(corr)<-shift;
  mn<-apply(corr, 2, function(r) weighted.mean(r, w));
  peak.ind<-which(mn>c(mn[-1], 1)&mn>c(1, mn[-length(mn)])&mn>quantile(mn)[4]);
  strand.correlation<-list(shift=shift, coefficient=mn, peak.index=peak.ind, all=corr);
  ############################################################################################################################### 
  
  ################################################## summarize peaks ################################################## 
  if (verbose) print('Summarizing ChIP-seq peaks ...');
  if (!stats$mapping$pair$paired) {
    if (frag.len<=0) frag.len<-shift[peak.ind[min(length(peak.ind), 2)]]; # if fragment length not given, make the best guess
    suppressWarnings(cov<-lapply(reads, function(r) coverage(resize(r, frag.len))[[1]]));
  } else {
    # TODO: paired end reads
  }
  
  # identify islands
  H<-ceiling(6*stats$depth$average);
  islands<-lapply(cov, function(r) slice(r, lower=H));
  islands<-lapply(names(reads), function(nm) GRanges(seqnames=nm, IRanges(start(islands[[nm]]), end(islands[[nm]])), height=max(islands[[nm]]), seqlengths=seqlengths(reads[[nm]])));
  suppressWarnings(peaks<-do.call('c', islands));
  sz<-genome$size[[1]];
  names(sz)<-rownames(genome$size);
  seqlevels(peaks)<-MapChrNames(seqlengths(peaks), sz);
  
  # Height summary
  hgt<-elementMetadata(peaks)[[1]];
  wth<-width(peaks);
  height<-c(sort(height[height>=H&height<max(hgt)]), max(hgt));
  ct<-sapply(height, function(h) length(hgt[hgt>=h]));
  mean.wd<-sapply(height, function(h) mean(wth[hgt>=h]));
  peak<-list(summary=data.frame("Height"=height, Count=ct, 'Average_width'=mean.wd));
  
  # Height distribution
  c<-c(1:100, seq(110, 200, 10), seq(300, 1000, 100), seq(2000, 10000, 1000));
  c<-c(c[c>=H&c<max(hgt)], max(hgt));
  n<-sapply(c, function(c) length(hgt[hgt>=c]));
  peak$height<-list(cutoff=c[n>0], count=n[n>0]);
  
  # Peak width
  dens<-density(log10(wth));
  dens$x<-exp(dens$x*log(10));
  peak$width<-list(summary=summary(wth), density=dens);

  c<-sapply(genome$features, function(f) countOverlaps(peaks, f, ignore.strand=TRUE)); 
  c[c>1]<-1;
  ct<-sapply(height, function(h) apply(c, 2, function(c) sum(c[hgt>=h])));
  colnames(ct)<-paste('Height>=', height, sep='');
  ct<-ct[, colSums(ct)>0];
  peak$count.by.feature<-ct;
  peak$peaks.per.mb<-10^6*ct/genome$total.feature.length;
  peak$all<-peaks;
  ############################################################################################################################### 
  
  ################################################## regions around TSS ######################################################### 
  if (verbose) print('Summarizing TSSs ...');
  dep<-lapply(reads, function(r) SummarizeTss(r, genome$tss, resize=1, max.dep=1));
  ss<-do.call('rbind', lapply(dep, function(x) x[[1]])); # sense strand
  as<-do.call('rbind', lapply(dep, function(x) x[[2]])); # antisense strand
  
  avg<-nrow(ss)*sum(as.numeric(bamchop$stats$depth$by.chromosome$Unique_mapping))/sum(as.numeric(bamchop$stats$depth$by.chromosome$Effecitive_size));
  
  tss<-list(global.average=avg, by.location=data.frame(Sense=colSums(ss), Antisense=colSums(as)), 
            by.gene=data.frame(RefSeq_ID=rownames(ss), Sense=rowSums(ss), Antisense=rowSums(as), Total=rowSums(ss)+rowSums(as)));
  ############################################################################################################################### 
  
  param<-list('fragment.length'=frag.len);
  
  bamchop$chip<-list(strand.correlation=strand.correlation, peak=peak, tss=tss, parameters=param);
  bamchop;
}