GetRegionsStats<-function(bam.file, chr, start, end, strand, names) {
  # bam.file              String, path to bam.file
  # chr,start,end,strand  query regions within which reads mapped to will be retrieved from bam file
  # read.len              maximum read length
  library(Rsamtools);
  library(R.oo);
  
  # use chromosome names as in bam file
  x<-chr %in% names(scanBamHeader(bam.file)[[1]][[1]]); 
  regions<-GRanges(seqnames=chr[x], ranges=IRanges(start[x], end[x]), strand=strand[x]);
  names(regions)<-names[x];
  seqlengths(regions)<-scanBamHeader(bam.file)[[1]][[1]][seqlevels(regions)];

  stat<-list();
  stat$index<-(1:length(x))[x]; # indexes of regions 
  
  if (sum(as.numeric(width(regions)))>=10^8) scan('GetReadStats::Warning: the query regions have a total size over 100 million bp; the scanBam process will probably take a long time');
  if (length(regions)>=500000) scan('GetReadStats::Warning: there are more than 500 thousand query regions; the scanBam process will probably take a long time');

  # read in 1000 regions at a time
  w<-mean(width(regions));
  if (w>10000) in.one<-5000000 else in.one<-2000000;
  l<-cumsum(width(regions));
  c<-ceiling(l[length(l)]/in.one);
  ind1<-sapply(1:c, function(c) (which(l>=(c*in.one)))[1]);
  ind1[length(ind1)]<-length(regions);
  ind0<-c(1, ind1[-length(ind1)]+1);
  
  counts<-list();
  mapqs<-list();
  bases<-list();
  quals<-list();

  for (i in 1:length(ind0)) {
    #print(i*in.one);
    #-----------------------------------------------------------------------------------------------#
    reads<-scanBam(bam.file, param=ScanBamParam(which=regions[ind0[i]:ind1[i]], what=scanBamWhat()));
    #-----------------------------------------------------------------------------------------------#
    names(reads)<-names[ind0[i]:ind1[i]];
    
    counts[[i]]<-t(sapply(reads, function(x) table(x[['strand']])));
    bases[[i]]<-t(sapply(reads, function(x) colSums(alphabetFrequency(x[['seq']]))));
    mapqs[[i]]<-lapply(reads, function(x) x[['mapq']]);
    
    qual<-lapply(reads, function(x) as.character(x[['qual']]));
    n<-sapply(qual, length);
    qual<-unlist(qual, use.names=F);
    len<-max(nchar(qual));
    
    sm<-rep(0, length(qual));
    ct<-rep(0, length(qual));
    
    for (j in 1:len) {
      sc<-charToInt(substring(qual, j, j))-33;
      sm[sc>=0]<-sm[sc>=0]+sc[sc>=0];
      ct[sc>=0]<-ct[sc>=0]+1;
    }
    
    quals[[i]]<-split(sm/ct, rep(1:length(reads), n));
    names(quals[[i]])<-names(reads)[n>0];
}

  stat<-list();
  stat$count<-do.call('rbind', counts);
  stat$mapq<-do.call('c', mapqs);
  stat$base<-do.call('rbind', bases);
  
  q<-do.call('c', quals);
  q[setdiff(names(regions), names(q))]<-NA;
  stat$qual<-q[names(regions)];
    
  stat$summary<-list(
    count=colSums(stat$count)[c('+', '-')],
    base=colSums(stat$base)[c('A', 'C', 'G', 'T', 'N')],
    mean.mapq=sapply(stat$mapq, function(x) mean(x, na.rm=TRUE)),
    mean.qual=sapply(stat$qual, function(x) mean(x, na.rm=TRUE))
  )
  
  stat;
}
