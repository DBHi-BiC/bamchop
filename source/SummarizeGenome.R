# Retrieve QC matrices from whole genome
SummarizeGenome <-function(bamfile, genome, include.chr=10^7, subset.size=10^5, verbose=TRUE) {
# bamfile       full path to the location of BAM file; named by sample name
# genome        Generic genome statistics, save in the /databases subfolder
# include.chr   names of chromosomes to be includesd in the report if string vector; minimum size of chromosomes to be reported if integer
# subset.size   number of randomly selected reads to be used to generate some of the QC matrices, such as sequencing quality. Raising this value linearly increases running time
  if (verbose) print('Summarizing genome ...');
  
  library(Rsamtools);
  library(R.oo);
  
  
  # check if the index file exists
  baifile<-paste(bamfile, '.bai', sep='');
  if (!file.exists(baifile)) {
    warning('the index file corresponding to the Bam file "', bamfile, '" does not exist. Bamchop is trying to generate one, which may take a long time');
    indexBam(bamfile);
  }
  
  # TODO: summarize SAM header
  
  chr<-scanBamHeader(bamfile)[[1]]$targets;
  if (class(include.chr)=='character')  chr<-chr[include.chr] else chr<-chr[chr>=include.chr];
  chr<-chr[!is.na(chr)];
  if (length(chr)==0) stop('Error: no valid chromosome to be included in the report');

  ##################################
  # read in mapping location of reads
  if (verbose) print('Reading in mapping locations ...');
  gr<-lapply(names(chr), function(x) GRanges(x, IRanges(1, chr[x]), seqlengths=chr[x]));
  names(gr)<-names(chr);
  #reads<-lapply(gr, function(gr) scanBam(bamfile, param=ScanBamParam(what=c('pos', 'qwidth', 'strand'), which=gr, flag=scanBamFlag(isUnmappedQuery=FALSE)))[[1]]);
  reads<-lapply(gr, function(gr) readBamGappedAlignments(bamfile, param=ScanBamParam(which=gr, flag=scanBamFlag(isUnmappedQuery=FALSE, isNotPrimaryRead=FALSE))));
  reads<-reads[sapply(reads, length)>0];
  chr<-chr[names(chr) %in% names(reads)];
  gapped<-sapply(reads, function(r) r[ngap(r)>0]);  # useful to RNAseq data
  cigar<-unlist(lapply(reads, cigar), use.names=FALSE);
  qwidth<-unlist(lapply(reads, qwidth), use.names=FALSE);
  reads<-lapply(names(chr), function(c) GRanges(c, IRanges(start(reads[[c]]), end(reads[[c]])), strand=strand(reads[[c]]), seqlengths=chr[c])[ngap(reads[[c]])==0])
  
  max.len<-max(sapply(reads, function(x) max(width(x))));
  cov<-lapply(reads, function(x) coverage(x)[[1]]);
  names(cov)<-names(reads)<-names(chr);
  
  # retrieve "flag" from BAM file, no filtering
  if (verbose) print('Reading in mapping flag ...');
  sam.flag<-unlist(lapply(gr, function(gr) scanBam(bamfile, param=ScanBamParam(what=c('flag'), which=gr))[[1]][[1]]), use.names=FALSE);
  
  #bamtools.stats<-BamtoolsStats(bamfile, paste(bamchop.path, '/bamtools/bamtools', sep=''));
  
  # number of reads mapped to each chromosome; three columns in order: chromosome length; mapped reads; and unmapped reads
  #samtools.stats<-SamtoolsStats(bamfile, paste(bamchop.path, '/samtools/samtools', sep=''))[names(chr), ];
  #total.mapped<-sum(samtools.stats[,2]);
  #total.unmapped<-sum(samtools.stats[,3]);
  
  # randomly select a subset of reads to generate QC stats using bamtools
  if (verbose) print('Reading in random read subset ...');
  ct.chr<-sapply(reads, length); 
  ns<-pmax(1, ceiling(ct.chr*subset.size/sum(as.numeric(ct.chr))));
  names(ns)<-names(chr); 
  selected<-lapply(names(chr), function(chr) if (ns[chr]<ct.chr[chr]) SampleBam(bamfile, cov[[chr]], ns[chr], chr) 
                   else scanBam(bamfile, param=ScanBamParam(what=scanBamWhat(), flag=scanBamFlag(isUnmappedQuery=FALSE)))[[1]]);
#  for (i in 1:length(selected)) {
#	selected[[i]][['seq']]<-as.character(selected[[i]][['seq']]);
#        selected[[i]][['qual']]<-as.character(selected[[i]][['qual']]);
#  }
  selected<-lapply(scanBamWhat(), function(fld) do.call('c', lapply(selected, function(x) as.vector(x[[fld]]))));
  names(selected)<-scanBamWhat();

  #remove duplicated entries
  qname<-selected[['qname']];
  rdm<-sample(1:length(qname), length(qname));
  dup<-duplicated(qname[rdm]);
  selected<-lapply(selected, function(x) x[rdm][!dup]);
  
  ############################## DONE with retrieving data from BAM file #################################################################
  
  # Summarizi ng ... 
  print('Summarizing statistics ...');
  stats<-list();
  stats$depth<-SummarizeDepth(cov, reads, genome);
  stats$quality<-SummarizeQuality(selected[['qual']], selected[['strand']]);
  stats$mapping<-SummarizeMapping(reads, sam.flag, qwidth, selected[['mapq']], cigar, selected[['isize']]);
  stats$bases<-SummarizeBases(selected[['seq']], selected[['strand']], genome);
  
  stats$summary<-c(
    'Number of chromosomes' = length(chr),
    'Total reference size (bp)' = sum(as.numeric(chr)),
    'Total effective size (bp)' = sum(as.numeric(genome$size$Effective)),
    'Total entries' = stats$mapping$total,
    'Total mapped reads' = stats$mapping$total.mapped,
    'Total upmapped reads' = stats$mapping$total.unmapped,
    'Total mappings' = stats$mapping$total-stats$mapping$total.unmapped,    
    'Total mapping locations' = sum(stats$mapping$duplicates$Location_count),
    
    'Base N%' = stats$bases$N[1, 'Percentage'],
    '(G+C)%' = stats$bases$observed.vs.expected[2, 'GC'], 
    'Mapped to forward strand%' = 100-stats$mapping$flag$by.bit['0X10', 2],
    'Duplicated mapping reads%' = 100-stats$mapping$duplicates[stats$mapping$duplicates[,1]==1, 4],
    
    'Best sequencing quality' = as.numeric(stats$qual$summary[6]),
    'Average sequencing quality' = as.numeric(stats$qual$summary[4]),

    'Maximum mapping length (bp)' = as.numeric(stats$mapping$length[6]),
    'Minimum mapping length (bp)' = as.numeric(stats$mapping$length[1]),   
    'Average mapping length (bp)' = as.numeric(stats$mapping$length[4]),

    'Best mapping quality' = as.numeric(stats$mapping$mapq$summary[6]),
    'Average mapping quality' = as.numeric(stats$mapping$mapq$summary[4]), 
    
    'Highest sequencing depth' = stats$depth$maximum, 
    'Average sequencing depth' = stats$depth$average,
    'Mapped reads per kilobase' = 1000*stats$mapping$total.mapped/sum(as.numeric(genome$size$Effective))
  )
  
  list(
    sample.name=names(bamfile),
    genome=genome,
    bamfile=file.info(bamfile),
    mapping.location=reads,
    coverage=cov,
    random.subset=selected,
    stats=stats
  );
}
