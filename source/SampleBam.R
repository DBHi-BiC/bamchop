# Randomly selected given number of reads from BAM file 
SampleBam<-function(bamfile, cov, n, chr.name) {
  n<-max(1, n);
  
  summ<-quantile(cov[cov>0]);
  up<-summ[4]+3*(summ[4]-summ[2]); # the maximum number of reads can come from one loc

  # number of regions to retrieve reads
  ind<-which(cov>0&cov<=up);
  n0<-ceiling(n/summ[3]);

  if (length(n0)>=length(ind)) loc<-ind else loc<-sample(ind, n0);
  g<-GRanges(chr.name, IRanges(loc, loc));  
  rd<-scanBam(bamfile, param=ScanBamParam(what=scanBamWhat(), which=g, flag=scanBamFlag(isUnmappedQuery=FALSE)));
#  for (i in 1:length(rd)) {
#    rd[[i]][['seq']]<-as.character(rd[[i]][['seq']]);
#    rd[[i]][['qual']]<-as.character(rd[[i]][['qual']]);
#  }
  rd<-lapply(scanBamWhat(), function(fld) unlist(lapply(rd, function(x) {
	if (fld=='seq' | fld=='qual') as.character(x[[fld]])
	else as.vector(x[[fld]]);
  }), use.names=FALSE));
  names(rd)<-scanBamWhat();
  #dup<-duplicated(rd[['qname']]);
  #rd<-lapply(rd, function(a) a[!dup]);

  rd;
}
