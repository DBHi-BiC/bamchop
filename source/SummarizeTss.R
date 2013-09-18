# summarize depth around TSS
SummarizeTss<-function(reads, tss, resize=0, max.dep=0, by.strand=TRUE, ws=1000) {
# reads   GRanges, read location on a chromosome
# tss     GRanges, TSS location
  ws<-max(0, ws);
  
  q<-seqlengths(reads)[as.vector(seqnames(reads[1]))];
  chr<-MapChrNames(q, seqlengths(tss)); 
  if (is.na(chr)) NA else { 
    tss<-tss[seqnames(tss)==chr]; 
    if (length(tss)>0) tss<-resize(tss, 1);
    if (length(tss)>0) tss<-tss[(start(tss)-ws)>0&(start(tss)+ws)<=q]; 
    if (length(tss)==0) NA else {      
      if (resize>0) reads<-resize(reads, resize);
      if (by.strand) cov<-list(forward=coverage(reads[strand(reads)=='+'])[[1]], reverse=coverage(reads[strand(reads)=='-'])[[1]]) else cov<-list(coverage=coverage(reads)[[1]]);
    
      str<-as.vector(strand(tss));
      loc<-start(tss);
      ref<-as.vector(elementMetadata(tss)$tx_name);
      
      d<-lapply(cov, function(cov) seqselect(cov, start=loc-ws, end=loc+ws));
      d<-lapply(d, as.numeric);
      d<-lapply(d, function(d) matrix(d, byrow=TRUE, nc=1+2*ws, nr=length(loc)));
      
      if (by.strand) {
        str<-as.vector(strand(tss));
        ss<-d[[1]];
        ss[str=='-']<-d[[2]][str=='-', ncol(d[[2]]):1];
        as<-d[[2]];
        as[str=='-']<-d[[1]][str=='-', ncol(d[[1]]):1];
        d<-list(sense=ss, antisense=as);
      }
      
      for (i in 1:length(d)) {
        if (max.dep>0) d[[i]][d[[i]]>max.dep]<-max.dep;
        rownames(d[[i]])<-ref;
        colnames(d[[i]])<-(-1*ws):ws;
      }
      
      d;
    }
  }
}