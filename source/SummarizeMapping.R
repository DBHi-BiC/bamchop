# Summarize mapping information, including SAM fields: flag, mapq, and cigar
SummarizeMapping<-function(reads, flag=NA, qwidth=NA, mapq=NA, cigar=NA, isize=NA, cutoffs=c(1:5, 10, 20, 30, 40)) {
# flag, mapq, cigar   SAM fields related to read mapping
# cutoffs mapq cutoffs

  mapping<-list(length=summary(qwidth), length.count=table(qwidth), total=length(flag));
  
  flag<-flag[!is.na(flag)];
  mapq<-mapq[!is.na(mapq)];
  cigar<-cigar[!is.na(cigar)];
  isize<-isize[!is.na(isize)];
  
  #################### SUMMARIZE MAPPING LENGTH ################################################
  #if (!identical(NA, seq) & length(seq)>0) {
  #  l<-nchar(seq);
  #  mapping$length<-list(total=length(l), all=l, summary=summary(l), density=density(l));
  #}
  
  #################### SUMMARIZE FLAG ################################################
  if (!identical(NA, flag) & length(flag)>0) {
    
      # defined by SAM format specifications
      desc<-c(
          '0X1'='template having multiple segments in sequencing',
          '0X2'='each segment properly aligned according to the aligner',
          '0X4'='segment unmapped',
          '0X8'='next segment in the template unmapped',
          '0X10'='SEQ being reverse complemented',
          '0X20'='SEQ of the next segment in the template being reversed',
          '0X40'='the first segment in the template',
          '0X80'='the last segment in the template',
          '0X100'='secondary alignment',
          '0X200'='not passing quality controls',
          '0X400'='PCR or optical duplicate');

  ct<-table(flag);

  smb<-c('00'='-', '01'='X'); # symbols used to represent yes (X) or no (-)
  
  # bits that represent the decimal values of the flag
  bits<-t(sapply(as.numeric(names(ct)), function(i) smb[as.character(intToBits(i)[1:length(desc)])]));
  bits<-data.frame(bits);
  rownames(bits)<-names(ct); 
  colnames(bits)<-names(desc);
  
  # summarize stats
  summ<-data.frame('Count'=as.vector(ct), 'Percentage'=as.vector(100*ct/length(flag)));
  summ<-cbind(summ, bits); 
  
  by.bit<-data.frame(Count=sapply(colnames(bits), function(b) sum(summ[[1]][bits[,b]=='X'])),
    Percentage=sapply(colnames(bits), function(b) sum(summ[[2]][bits[,b]=='X'])))
  
  mapping$total.unmapped<-by.bit['0X4', 'Count']; 
  mapping$total.mapped<-sum(as.numeric(sapply(reads, length))); 
  mapping$total.forward<-100-100*by.bit['0X10', 'Percentage'];  
      
  mapping$flag<-list(total=length(flag), description=data.frame(Description=desc), summary=summ, by.bit=by.bit);
  }

  #################### SUMMARIZE MAPQ ################################################
  if (!identical(NA, mapq) & length(mapq)>0) {
    summ<-summary(mapq);
    dist<-table(mapq);
    dist<-data.frame(Score=as.numeric(names(dist)), Count=as.numeric(dist));
    
    # number of reads with mapping score higher that given cutoffs
    cutoffs<-cutoffs[cutoffs>0&cutoffs<summ[6]];
    if (length(cutoffs)>0) ct<-c(sapply(cutoffs, function(c) sum(as.numeric(dist[dist[,1]>=c,2]))), dist[nrow(dist),2]) else ct<-dist[nrow(dist), 2];
    if (dist[1,1]==0) ct<-c(dist[1,2], ct) else ct<-c(0, ct);
    pct<-round(100*ct/sum(as.numeric(dist[,2])), 2);
    count<-data.frame(Count=ct, Percentage=pct);
    rownames(count)<-c('mapq=0', paste('mapq>=', cutoffs, sep=''), paste('mapq=', summ[6], sep=''));
    
    mapping$mapq<-list(total=length(mapq), all=mapq, summary=summ, distribution=dist, count=count);
  }
  
  #################### SUMMARIZE CIGAR ################################################
  if (!identical(NA, mapq) & length(cigar)>0) {
    # defined by SAM format specifications
    desc<-c(
      "M"="alignment match (can be a sequence match or mismatch)",
      "I"="insertion to the reference",
      "D"="deletion from the reference",
      "N"="skipped region from the reference",
      "S"="soft clipping (clipped sequences present in SEQ)",
      "H"="hard clipping (clipped sequences NOT present in SEQ)", 
      "P"="padding (silent deletion from padded reference)",
      "="="sequence match",
      "X"="sequence mismatch"
    );
    
    ct<-sapply(names(desc), function(x) length(cigar[grep(x, cigar)]));
    ct<-ct[ct>0];
    pct<-round(100*ct/length(cigar), 4);
    count<-data.frame(Count=ct, Percentage=pct);
    
    gp<-cigar[grep('N', cigar)];
    if (length(gp)>0) {
      len<-unlist(lapply(gp, function(x) ParseCigar(x, 'N')), use.names=FALSE);
      dens<-density(log10(len));
      dens$x<-exp(dens$x*log(10));
      gapped<-list(count=length(gp), summary=summary(len), density=dens);
    }
    else gapped<-list(count=0);
    
    mapping$cigar<-list(total=length(cigar), description=data.frame(Description=desc), count=count, gapped=gapped);
  }
  
  #################### SUMMARIZE DUPLICATION ################################################  
  max.reads<-10^7
  ct<-lapply(reads, function(gr) {
    gr<-resize(gr, 1);
    loc<-start(gr);
    str<-as.vector(strand(gr))=='-';
    loc[str]<--1*loc[str];
    as.numeric(table(loc));
#    str<-as.vector(strand(gr));
#    fw<-start(gr)[str=='+'|str=='*']; # reads mapped to forward or unknown strand
#    as.numeric(table(fw[1:min(max.reads, length(fw))]));
  })
  
  summ<-table(unlist(ct, use.names=FALSE));
  lvl<-as.numeric(names(summ));
  ct<-lvl*as.numeric(summ);
  mapping$duplicates<-data.frame(Level=lvl, Location_count=as.numeric(summ), Read_count=ct, Percentage=as.numeric(100*ct/sum(ct)));
  
  #################### SUMMARIZE PAIRS ################################################
  pair<-list(paired=by.bit['0X1', 1]>0);
  if(pair$paired) {
    pair$summary<-by.bit[c('0X1', '0X2', '0X4'), ];
    rownames(pair$summary)[1:3]<-c('Total paired-end reads', 'Ends properly mapped', 'One end unmapped');
    isize<-abs(isize[!is.na(isize) & isize!=0]);
    if (length(isize)>1) {
      summ<-summary(isize);
      dens<-density(log10(isize));
      dens$x<-exp(dens$x*log(10));
      pair$insert.size<-list(summary=summ, density=dens);
    }
  }
  mapping$pair<-pair;
  
  mapping;
}