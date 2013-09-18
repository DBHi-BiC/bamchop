# count number of reads mapped to unique locations
CountNoDup<-function(reads) {
  r1<-reads[strand(reads)=='+'];
  r2<-reads[strand(reads)=='-'];
  
  l1<-start(r1);
  l2<-end(r2);
  
  length(unique(l1))+length(unique(l2));
}

# Get position-specfic average depth of a set of regions
AvgDepth<-function(cov, n.int=100) {
  # cov     a list of Rle, each is the depth data within a genomic regions
  # n.int   number of intervals to split each region
  
  n<-sapply(cov, length);
  if(n.int<=0 | n.int>max(n)) n.int<-max(n); 
  
  dep<-unlist(lapply(cov, function(x) as.vector(x)/mean(x)), use.names=FALSE);
  loc<-round(n.int*unlist(lapply(n, function(n) ((1:n)-0.5)/n), use.names=FALSE));
  wgt<-rep(1/n, n)

  x<-split(dep, loc);
  y<-split(wgt, loc);
  
  z<-rep(NA, n.int+1);
  names(z)<-0:n.int;
  
  for (i in 1:length(x)) z[[names(x)[i]]]<-weighted.mean(x[[i]], y[[i]], na.rm=TRUE); 
  
  z;
}

# retrieve depth information of given segments from all chromosomes
GetAllSegments<-function(values, ranges) {
  # values      A list of numeric vectors or Rles
  # ranges      an GRanges objects defining the ranges of data to be selected
  
  chr<-seqlevels(ranges);
  chr<-chr[chr %in% names(values)];
  if (length(chr)<1) stop('No valid chromosome names.');
  
  # index of segments on each chromosome
  ind<-lapply(chr, function(chr) which(seqnames(ranges)==chr));
  names(ind)<-chr;
  ind<-ind[sapply(ind, length)>0];
  chr<-names(ind);
  
  segments<-lapply(chr, function(chr) {
    g<-ranges[ind[[chr]]];
    segments<-GetSegments(values[[chr]], start(g), end(g));
    names(segments)<-ind[[chr]];
    segments;
  });
  
  segments<-do.call('c', segments);
  segments[order(as.numeric(names(segments)))];
}


# Get segments of non-fixed length on a long vector (ex. depth on a chromosome) around a given set of indexes
# return a list of numeric vectors  with each element represents a segment
GetSegments<-function(v, start, end, as.vector=FALSE) {
  # v; a numerical vector or Rle object from which the segments are retrieved
  # start; Integer vector, start locations
  # end; Integer vector, end locations, same length as 'start' 
  # as.vector; if TRUE, convert Rle elements to vectors
#  library(chipseq);
  
  start[start<1]<-1; start[start>length(v)]<-length(v); 
  end[end<1]<-1; end[end>length(v)]<-length(v); 
  
  x<-seqselect(v, start=start, end=end); # retrieve subset into a single vector
  if (as.vector) x<-as.vector(x); 
  
  size<-end-start+1; # size of ranges
  factor<-rep(1:length(size), size); # splitting factor
  
  segment<-split(x, factor);
  segment;
}
