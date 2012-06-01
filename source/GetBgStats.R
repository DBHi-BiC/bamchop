# Summarize statistics of selected regions of background genome
GetBgStats<-function(bam.file, chr, start, end, strand) {
# bam.file              String, path to bam.file
# chr,start,end,strand  query regions within which reads mapped to will be retrieved from bam file
library(Rsamtools);
library(R.oo);

# use chromosome names as in bam file
x<-chr %in% names(scanBamHeader(bam.file)[[1]][[1]]); 
regions<-GRanges(seqnames=chr[x], ranges=IRanges(start[x], end[x]), strand=strand[x]);

seqlengths(regions)<-scanBamHeader(bam.file)[[1]][[1]][seqlevels(regions)];
names(regions)<-(1:length(x))[x];

stat<-list();
stat$index<-(1:length(x))[x]; # indexes of regions 

if (sum(as.numeric(width(regions)))>=10^8) scan('GetReadStats::Warning: the query regions have a total size over 100 million bp; the scanBam process will probably take a long time');
if (length(regions)>=500000) scan('GetReadStats::Warning: there are more than 500 thousand query regions; the scanBam process will probably take a long time');

#------------------------------------------------------------------------------#
reads<-scanBam(bam.file, param=ScanBamParam(which=regions, what=scanBamWhat()));
#------------------------------------------------------------------------------#
n<-sapply(reads, function(x) length(x[[1]])); 
values<-list(
  id=rep(names(regions), n),
  chr=unlist(lapply(reads, function(x) x[['rname']]), use.names=FALSE),
  pos=unlist(lapply(reads, function(x) x[['pos']]), use.names=FALSE),
  str=unlist(lapply(reads, function(x) as.vector(x[['strand']])), use.names=FALSE),
  mapq=unlist(lapply(reads, function(x) x[['mapq']]), use.names=FALSE),
  cigar=unlist(lapply(reads, function(x) x[['cigar']]), use.names=FALSE),
  seq=unlist(lapply(reads, function(x) as.character(x[['seq']])), use.names=FALSE),
  qual=unlist(lapply(reads, function(x) as.character(x[['qual']])), use.names=FALSE)
); 

# in case unmapped reads are included, exclude those reads
pos<-values[['pos']]; 
values<-lapply(values, function(v) v[!is.na(pos)]);

# split into forward vs. backward reads
str<-values[['str']];
fw<-lapply(values, function(v) v[str=='+']);
bw<-lapply(values, function(v) v[str=='-']);


# Total counts
stat$counts<-c(total=length(str), forward=length(fw[[1]]), backward=length(bw[[1]]));

# Distribution of mapping scores
mapq<-table(values[['mapq']]);
mapq<-mapq[mapq>0];
stat$mapq<-mapq[order(as.numeric(names(mapq)))];

# Information from CIGAR string
cigar<-values[['cigar']];
cigar<-sapply(c('I', 'D', 'S', 'H'), function(c) length(grep(c, cigar)));
names(cigar)<-c('Insertion', 'Deletion', 'Soft clipping', 'Hard clipping');
stat$cigar<-cigar[cigar>0];

# position specific base frequecy
bases<-list();
seq<-c(fw[['seq']], as.character(reverseComplement(DNAStringSet(bw[['seq']])))); # reverse backward reads
len<-max(nchar(seq)); # maximum length
freq<-matrix(nr=len, nc=5);
colnames(freq)<-c('A', 'C', 'G', 'T', 'N');
rownames(freq)<-1:len;
for (i in 1:len) freq[i, ]<-table(substr(seq, i, i))[colnames(freq)];
freq[is.na(freq)]<-0;
bases$frequency<-freq;

# counts of combinations of first 2 bases, Ns excluded
di<-table(substr(seq, 1, 2));
bases$di.base<-di[-grep('N', names(di))];

stat$bases<-bases;

# position specific quality score
# TODO: also use backward reads
qual<-fw[['qual']];
score<-lapply(1:len, function(i) substr(qual, i, i));
score<-lapply(score, function(x) x[x!='']);
score<-lapply(score, function(x) charToInt(x)-33); # convert to scores, use the convention of Sanger sequencing
names(score)<-1:len;
mean<-sapply(score, function(x) if(length(x)==0) NA else mean(x));
sd<-sapply(score, function(x) if (length(x)<2) 0 else sd(x)); 
stat$quality$summary<-cbind(position=1:len, mean=mean, sd=sd, n=sapply(score, length));
stat$quality$density<-density(unlist(score, use.names=FALSE), bw=1);
stat$quality$percentiles<-sapply(score, function(s) quantile(s, probs=seq(0, 1, 0.01))); 

stat;
}
