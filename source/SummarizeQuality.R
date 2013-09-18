# summarize sequencing quality
SummarizeQuality<-function(q, strand=NA, adjust=-33, cutoffs=c(5, 10, 13, 20, 30, 40)) {
# q       Vector of string providing position specific quality score of reads
# strand  Strand the reads was mapped to; must have the same length as q
# adj     The adjustment of the score after converting the ASCII code to Phred score; the default is -33 if the convention of Sanger sequencing is used by the aligner
  
  if (class(q) != 'character') q<-as.character(q);

  L<-nchar(q); # length

  # score of individual reads
  q.read<-lapply(strsplit(q, ''), function(q) charToInt(q)+adjust);
  by.read<-data.frame(Length=L, 'Mean_score'=sapply(q.read, function(q) mean(q, na.rm=TRUE)));

  # all scores
  scores<-unlist(q.read, use.names=FALSE);
  
  # reverse order of reads mapped to - strand
  if (!identical(strand, NA) & length(strand)==length(q)) q.read[strand=='-']<-lapply(q.read[strand=='-'], function(x) x[length(x):1]);
  
  # position specific scores
  len<-max(L);
  score<-lapply(1:len, function(i) sapply(q.read[L>=i], function(x) x[i]));
  score<-lapply(score, function(x) x[!is.na(x)]); 
  names(score)<-1:len;
  mean<-sapply(score, function(x) if(length(x)==0) NA else mean(x));
  n<-sapply(score, length);
  sd<-sapply(score, function(x) if (length(x)<2) 0 else sd(x)); 
  
  # Score counts by cutoff values
  Max<-max(scores, na.rm=TRUE);
  cutoffs<-sort(cutoffs[cutoffs>0&cutoffs<Max]);
  if (length(cutoffs)>0) ct<-sapply(cutoffs, function(c) length(scores[scores>=c])) else c();
  ct<-c(length(scores[scores==0]), ct, length(scores[scores==Max]));
  pct<-round(100*ct/length(scores), 2);
  count<-data.frame(Count=ct, Percentage=pct);
  rownames(count)<-c('Score=0', paste('Score>=', cutoffs, sep=''), paste('Score=', Max, sep=''));
  
  list(
  total.reads=length(q),
  total.bases=length(scores),
  all=scores,
  summary=summary(scores),  
  by.read=by.read,
  count=count,
  density=density(scores, bw=1),
  position.specific=list(summary=data.frame('Position'=1:len, 'Read_count'=n, 'Mean_score'=mean, 'Standard_deviation'=sd), percentiles=sapply(score, function(s) quantile(s, probs=seq(0, 1, 0.01))))
  )
}