BamtoolsRandom<-function(bamfile, bamtools.path, chr, n, verbose=FALSE) {
  library(R.utils); 
  
  if (n<=0) n<-1;
  seg<-ceiling(n/80); 
  files<-paste('random', 1:seg, '.bam', sep='');
  file.remove(files); 
  
  rd<-list();
  for (i in 1:length(files)) {
    if (verbose) cat('generating random reads ', i*100-99, ' to ', i*100, ' on chromosome ', chr, '\n');
    system(paste(bamtools.path, ' random -in ', bamfile, ' -out ', files[i], ' -region ', chr, ' -n 100', sep=''));
    rd[[i]]<-scanBam(files[i], param=ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isUnmappedQuery=FALSE)))[[1]];
  }
  file.remove(files);
  
  rd<-lapply(scanBamWhat(), function(fld) do.call('c', lapply(rd, function(x) x[[fld]])));
  names(rd)<-scanBamWhat();  
  dup<-duplicated(rd[['qname']]);  
  rd<-lapply(rd, function(a) a[!dup]);
  if(length(rd[[1]])>n) {
    ind<-sample(1:length(rd[[1]]), n);
    rd<-lapply(rd, function(a) a[ind]);
  }
  rd;  
}

# Get stats retrieved by bamtools
BamtoolsStats<-function(bamfile, bamtools.path) {
  system(paste(bamtools.path, ' stats -in ', bamfile, ' >stats.txt', sep=''));
  stats<-scan('stats.txt', what=list(''), sep='\t')[[1]][c(6, 8, 12, 14, 16, 20, 24)];
  stats<-sub(')', '', stats);
  stats<-sub('\\(', '', stats);
  names(stats)<-c('Mapped reads', 'Forward strand', 'Failed QC', 'Duplicates', 'Paired-end reads', 'Both ends mapped', 'Singletons')
  file.remove('stats.txt');
  stats;
}

# Get read count on each chromosome
SamtoolsStats<-function(bamfile, samtools.path) {
  system(paste(samtools.path, ' idxstats ', bamfile, ' >counts.txt', sep=''));
  count<-read.table('counts.txt', sep='\t', row=1);
  file.remove('counts.txt');
  count;
}

