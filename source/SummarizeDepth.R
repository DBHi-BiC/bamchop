# Summarizie depth information on multiple chromosomes
SummarizeDepth<-function(cov, reads, genome, cutoffs=c(1, 5, 10, 20, 30, 50, 100, 1000, 10000)) {
# cov		  A list of named Rle, Sequencing depth along each chromosome
# reads   Mapping location of reads
# genome	Genome metadata

# give temporay names if not available
names(cov)[is.na(names(cov))]<-(1:length(cov))[is.na(names(cov))];
  
chr.len<-sapply(cov, length);

# by genomic features such as exon, intron, ...
features<-genome$features;
f.dep<-sapply(features, function(f) {
  if (class(f) !='GRanges') NA else {
    cv<-coverage(f);
    chr.map<-MapChrNames(sapply(cov, length), sapply(cv, length));
    chr.map<-chr.map[!is.na(chr.map)];    
    if (length(chr.map)==0) NA else{
      c<-lapply(names(chr.map), function(x) cov[[x]][cv[[chr.map[x]]]>0]);
      weighted.mean(sapply(c, mean), sapply(c, length));
    }
  }
});


# remove assembly gaps before calculating stats
L<-sapply(cov, length); # full length of chromosomes
gaps<-genome$gaps; # assembly gaps of all chromosomes
if (class(gaps)=='GRanges'&length(gaps)>0) {
	cgaps<-coverage(gaps);
	lgaps<-sapply(cgaps, length); 
	cov<-lapply(names(cov), function(nm) {
		if (nm %in% names(cgaps)) c<-cgaps[[nm]] else c<-cgaps[lgaps==L[nm]][[1]];
		if (length(c)==0) cov[[nm]] else cov[[nm]][c==0]; 
	})
}

# depth by chromosomes
eff<-sapply(cov, length); # effective size
max<-sapply(cov, max);
avg<-round(sapply(cov, function(c) mean(as.numeric(c))),2);
max.loc<-lapply(cov, function(x) which(x==max(x)));
max.loc<-sapply(max.loc, function(x) x[ceiling(length(x)/2)]);
uni<-sapply(reads, CountNoDup);
by.chr<-data.frame('Chromosome_length'=L, 'Effecitive_size'=eff, 'Total_reads'=sapply(reads, length), 'Unique_mapping'=uni, 'Average_depth'=avg, 'Maximum_depth'=max, 'Maximum_location'=max.loc);

# Depth count
Max<-max(max);
cutoffs<-cutoffs[cutoffs>0&cutoffs<Max]; # make sure the first cutoff is 0
if (length(cutoffs)>0) ct<-sapply(cutoffs, function(c) sum(sapply(cov, function(cov) as.numeric(length(cov[cov>=c]))))) else c();
ct<-c(sum(sapply(cov,function(cov) as.numeric(length(cov[cov==0])))), ct, sum(sapply(cov,function(cov) as.numeric(length(cov[cov==Max])))));
pct<-round(100*ct/sum(as.numeric(eff)), 2);
count<-data.frame(Count=ct, Percentage=pct);
rownames(count)<-c('Depth=0', paste('Depth>=', cutoffs, sep=''), paste('Depth=', Max, sep=''));

# output object
list(average=weighted.mean(avg, w=eff), maximum=max(by.chr[['Maximum_depth']]), cutoffs=count, by.chromosome=by.chr, by.feature=f.dep);
}
