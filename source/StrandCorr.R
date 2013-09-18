# Calculate depth correlation between two strands
StrandCorr<-function(reads, resize=0, shift=0, nodup=TRUE) {
# reads   GRanges, location of mapped reads
# resize  Resize reads to given width; no resizing if 0
# shift   Non-negative integer vector, index shifting of the two strand
# nodup   If TRUE, remove duplicated mapping (max depth = 1)
  library(IRanges);
  
  if (resize>0) reads<-resize(reads, resize);
  fw<-coverage(reads[strand(reads)=='+'])[[1]];
  rv<-coverage(reads[strand(reads)=='-'])[[1]];
  
  if (nodup) {
    fw[fw>1]<-1;
    rv[rv>1]<-1;
  }
  
  list(weight=sum(fw)+sum(rv), coefficient=shiftApply(shift, fw, rv, cor));
}