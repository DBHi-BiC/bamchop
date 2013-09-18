# Summarize duplicated mapping
SummarizeDuplicates<-function(GRlist, max.reads=10^7) {
# GRlist    A list of GRanges objects, each represents mapping locations of reads of one chromosome; must be sorted by chromosomal locations
# max.reads If total number of reads is greater than this number on one chromosome, trim to this number to reduce time

  ct<-lapply(GRlist, function(gr) {
    str<-as.vector(strand(gr));
    fw<-start(gr)[str=='+'|str=='*']; # reads mapped to forward or unknown strand
    as.numeric(table(fw[1:min(max.reads, length(fw))]));
  })
  
  summ<-table(unlist(ct, use.names=FALSE));
  data.frame(Level=as.numeric(names(summ)), Count=as.numeric(summ), Percentage=as.numeric(100*summ/sum(summ)));
}