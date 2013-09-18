# Summarize nucleic acid bases in sequencing reads
SummarizeBases<-function(seq, strand=NA, genome=NA, end.n=10, k=5) {
# seq     Read sequences
# strand  Mapped to strand, must have the same length as seq  
# genome  Genome metadata
# end.n   The number of bases at both ends to summarize position-specific base frequency
# k       Check the frequency of the first and last k bases of all reads to evaluate possible bias and residual of adaptor comtamination
  
  if (class(seq) != "DNAStringSet") seq<-DNAStringSet(seq);
  
  if (!identical(strand, NA) & length(strand)==length(seq)) seq[strand=='-']<-reverseComplement(seq[strand=='-']);
  
  if (k>6) warning('The value of k=', k, ' for kmer checking is larger than 6, try to use smaller number to avoid prolonged computation.')
  if(end.n<k) end.n<-k;
  
  b4<-c('A', 'C', 'G', 'T');
  
  # overall frequency
  freq<-alphabetFrequency(seq);

  # summarize N base occurance
  N<-freq[, 'N'];
  summN<-data.frame(Total=c(sum(as.numeric(freq)), length(seq)), N=c(sum(N), length(N[N>0])));
  rownames(summN)<-c('Base', 'Read');
  summN$Percentage<-100*summN[,2]/summN[,1];
  
  # ACGT frequency of individual reads
  by.read<-freq[, b4];
  
  # Overall observed vs. expected frequency
  if (!identical(genome, NA)) {
    exp.ct<-colSums(genome$size[, colnames(by.read)]);
    exp<-c(exp.ct/sum(exp.ct), 'GC'=sum(exp.ct[c('G', 'C')])/sum(exp.ct));
  } else exp<-c(rep(1/length(b4), length(b4)), 0.5); 
  obs.ct<-colSums(by.read);
  obs<-c(obs.ct/sum(obs.ct), 'GC'=sum(obs.ct[c('G', 'C')])/sum(obs.ct));
  overall<-data.frame(100*rbind(exp, obs, obs/exp));
  rownames(overall)<-c('Expected(%)', 'Observed(%)', 'Observed/Expected(%)');
  
  # position-specific, head
  fst<-sapply(1:end.n, function(i) table(substr(seq, i, i))[b4]);
  fst[is.na(fst)]<-0;
  fst.ori<-fst/colSums(fst);
  fst.norm<-apply(fst, 2, function(x) 0.25*x/(obs[1:4]*sum(x)) )
  colnames(fst.norm)<-colnames(fst.ori)<-colnames(fst)<-1:end.n;
  
  # position-specific, tail
  rv<-reverse(seq);
  lst<-sapply(1:end.n, function(i) table(substr(rv, i, i))[b4]);
  lst[is.na(lst)]<-0;
  lst.ori<-lst/colSums(lst);
  lst.norm<-apply(lst, 2, function(x) 0.25*x/(obs[1:4]*sum(x)))[, ncol(lst):1];
  colnames(lst.norm)<-colnames(lst.ori)<-colnames(lst)<--1*(end.n:1);
  
  # counts of combinations of first 2 bases, Ns excluded
  di<-table(substr(seq, 1, 2))[paste(rep(b4, each=length(b4)), rep(b4, length(b4)), sep='')];
  di[is.na(di)]<-0;
  di.exp<-rep(fst.ori[,1], each=length(b4))*rep(fst.ori[,2], length(b4));
  names(di.exp)<-names(di);
  di.count<-data.frame(Count=di, "Expected_frequency"=di.exp, "Observed_frequency"=di/sum(di));
  di.frq<-t(matrix(di.count[,3]/di.count[,2], nr=4));
  dimnames(di.frq)<-list('First base'=c('A', 'C', 'G', 'T'), 'Second base'=c('A', 'C', 'G', 'T'));
  
  ################################ kmer frequency  ################################ 
  combs<-as.matrix(expand.grid(rep(list(1:length(b4)), k)));
  kmer<-apply(combs, 1, function(c) paste(b4[c], collapse=''));
  
  # first kmer at the beginning
  kmer.fst<-table(substr(seq, 1, k))[kmer];
  kmer.fst[is.na(kmer.fst)]<-0;
  fst.exp<-sum(kmer.fst)*apply(sapply(1:k, function(i) fst.ori[,i][combs[,i]]), 1, prod);
  fst.exp<-fst.exp*sum(kmer.fst)/sum(fst.exp);
  names(fst.exp)<-names(kmer.fst)<-kmer;
  k.first<-data.frame("Expected_count"=fst.exp, "Observed_count"=kmer.fst);
  k.first$"Observed/Expected"<-k.first[,2]/k.first[,1];
  rownames(k.first)<-paste('|', rownames(k.first), sep='');
  
  # last kmer at the end
  kmer.lst<-table(substr(rv, 1, k))[kmer];
  kmer.lst[is.na(kmer.lst)]<-0;
  names(kmer.lst)<-reverse(DNAStringSet(kmer));
  lst.exp<-sum(kmer.lst)*apply(sapply(1:k, function(i) lst.ori[,i][combs[,i]]), 1, prod);
  lst.exp<-lst.exp*sum(kmer.lst)/sum(lst.exp);
  names(lst.exp)<-names(kmer.lst)<-reverse(DNAStringSet(kmer));    
  k.last<-data.frame("Expected_count"=lst.exp, "Observed_count"=kmer.lst);
  k.last$"Observed/Expected"<-k.last[,2]/k.last[,1];
  rownames(k.last)<-paste(rownames(k.last), '|', sep='');
  
  k.outliers<-rbind(k.first[k.first$"Observed/Expected">3, ], k.last[k.last$"Observed/Expected">3, ]);
  k.outliers<-k.outliers[order(k.outliers$"Observed/Expected"), ][nrow(k.outliers):1, ];
  
  list(
    'N' = summN,
    'by.read' = by.read,
    "observed.vs.expected" = overall,
    "position.specific" = list(
      first=list(count=fst, freq.original=fst.ori, freq.normalized=fst.norm),
      last=list(count=lst, freq.original=lst.ori, freq.normalized=lst.norm),
      dibase=list(
        count=di.count,
        normalize.frequency=di.frq
      ),
      kmer=list(
        k=k,
        first=k.first,
        last=k.last,
        outliers=k.outliers
      )
    )
  );
}
