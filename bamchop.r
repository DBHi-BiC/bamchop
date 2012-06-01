library(Biostrings);
library(GenomicRanges);
library(Rsamtools);
library(seqLogo);

library(R.oo); # needed for ASCII
library(plotrix); # needed for pie3D() 
library(gplots);
library(vcd); # needed for mosaic()
library(xtable); # needed for xtable()

# load in custom functions
source('bamchop_param.r')
source('source/PlotLandscape.R');
source('source/PlotMisc.R');
source('source/Depth2BedGraph.R');
source('source/GetBgStats.R');
source('source/GetRegionsStats.R');


# bam file names
x<-scan(bams, flush=TRUE, sep='\t', what=list('', ''));
bam.list<-x[[2]];
names(bam.list)<-x[[1]];

# targeted regions
if (class(targets)=='character') {
  x<-scan(targets, flush=TRUE, sep='\t', what=list('', '', 0, 0, 0));
  targets<-GRanges(seqnames=x[[2]], ranges=IRanges(x[[3]], x[[4]]), toi=x[[5]]);
  names(targets)<-x[[1]];
} else if (class(targets)!='GRanges') {
stop("Class of variable <targets>", class(targets), "not recognized");
}
#targets<-eval(parse(text=load(targets.file)));
targets<-targets[order(as.vector(seqnames(targets)), start(targets), end(targets))];

if(!file.exists(path.out)) dir.create(path.out);

# Load reference genome, and select chromosomes to be included 
# By default, all autosomes, X, Y, and M will be included, random chromosome will not be included
g<-tolower(genome.name); 
if (g %in% c('hg19', 'human', 'homo sapiens', 'hs', 'h.s.', 'hsapiens', 'grch37')) {
  library(BSgenome.Hsapiens.UCSC.hg19);
  assign('ref.genome', Hsapiens);
  genome<-eval(parse(text=load('database/hg19.rdata')));
} else if (g %in% c('hg18', 'ncbi36')) {
  library(BSgenome.Hsapiens.UCSC.hg18);
  assign('ref.genome', Hsapiens);
  genome<-eval(parse(text=load('database/hg18.rdata'))); # TODO: create file
} else if (g %in% c('mm9', 'mouse', 'mouse musculus', 'mmusculus', 'ncbi37')) {
  library(BSgenome.Mmusculus.UCSC.mm9);
  assign('ref.genome', Mmusculus);
  genome<-eval(parse(text=load('database/hg19.rdata'))); # TODO: create file
} else if (g %in% c('mm8')) {
  library(BSgenome.Mmusculus.UCSC.mm8);
  assign('ref.genome', Mmusculus);
  genome<-eval(parse(text=load('database/hg18.rdata'))); # TODO: create file
} else stop ('Error: Unknown genome name', genome.name);
ref.chr<-genome$chr.names;
genome.size<-genome$genome.size;
gap.size<-genome$gap.size;
one.percent<-genome$one.percent[[1]];

# targets and reference genome have no common chromosome names
tg.chr<-seqlengths(targets);
if (length(intersect(names(tg.chr), names(ref.chr)))==0) {
  if (!map.chr.names) 
    stop("\tWe detect targets and reference genome share no common chromosome names.", 
         "\n\n\tWe recommend to re-format data using consistent names.", 
         "\n\n\tOtherwise, set value of <map.chr.names> to TRUE, and the program will use its best guess to map chromosome names.",
         "\n\n\tChromosome names of reference genome: ", paste(names(ref.chr), collapse='; '), 
         "\n\n\tChromosome names of targets: ", paste(names(tg), collapse='; '));
}

# if there is at least one chromosome name of targets not found in reference genome, map names and then remove names not in reference genome
if (min(seqlevels(targets) %in% names(ref.chr))==0) {
if (map.chr.names) {
  scan('!Warning: We noticed that there are chromosome names of the targeted regions not found in the reference genome. \nWe are trying to map the names using our guess. \nIn the future, we recommend to re-compile your target region file using names consistent to the reference.');
  map.tg.chr<-sapply(names(tg.chr), function(x) names(ref.chr)[ref.chr==tg.chr[x]][1]); 
} else {
  scan('!Warning: We noticed that there are chromosome names of the targeted regions not found in the reference genome. \nWe will only include those having common names from now on. \nIn the future, we recommend to re-compile your target region file using names consistent to the reference.');  
  map.tg.chr<-sapply(names(tg.chr), function(x) names(ref.chr)[names(ref.chr)==x][1]);    
}
targets<-targets[seqnames(targets) %in% names(map.tg.chr[!is.na(map.tg.chr)])]; 
nm<-names(targets);
targets<-GRanges(seqnames=map.tg.chr[as.vector(seqnames(targets))], range=IRanges(start(targets), end(targets)), strand=strand(targets)); 
names(targets)<-nm;
}
seqlengths(targets)<-ref.chr[seqlevels(targets)];

dates<-list();
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
for (i in 1:length(bam.list)) { dates[[i]]<-date();
  sample.name<-names(bam.list)[i]; print(sample.name);
  bam.file<-bam.list[i];
  if (verbose) cat('Processing sample', sample.name, '\n');
  
  # if there is inconsistent chromosome names
  bam.chr<-scanBamHeader(bam.file)[[1]][[1]];  
  if (length(intersect(names(bam.chr), names(ref.chr)))==0) {
    if (!map.chr.names) 
      stop("\tWe detect bam file and reference genome share no common chromosome names.", 
           "\n\n\tWe recommend to re-generate bam file using consistent names.", 
           "\n\n\tOtherwise, set value of <map.chr.names> to TRUE, and the program will use its best guess to map chromosome names.",
           "\n\n\tChromosome names of reference genome: ", paste(names(ref.chr), collapse='; '), 
           "\n\n\tChromosome names in bam file: ", paste(names(bam.chr), collapse='; '));
    }
  
  if (map.chr.names) {
    map.bam.chr<-sapply(names(bam.chr), function(x) names(ref.chr)[ref.chr==bam.chr[x]][1])
  } else {
    map.bam.chr<-sapply(names(bam.chr), function(x) names(ref.chr)[names(ref.chr)==x][1]);
  }
  map.bam.chr.rev<-names(map.bam.chr);
  names(map.bam.chr.rev)<-map.bam.chr;

  #@@@@@@
  #if (!exists('gr') | !exists('cov')) {
  cat('Reading data from', bam.file, '\n');
  #-------------------------------READ FROM BAM FILE-------------------------------------------#
  reads<-scanBam(bam.file, param=ScanBamParam(what=c('rname', 'strand', 'pos', 'qwidth')))[[1]];
  #--------------------------------------------------------------------------------------------#
  
  #----------------------------------------------------- Whole genome --------------------------------------------------------#
  # convert to GRanges object
  reads[['rname']]<-as.vector(map.bam.chr[reads[['rname']]]);
  vld<-!is.na(reads[['pos']]) & !is.na(reads[['rname']]); # whether the alignment is valid (has a specified location)
  gr<-GRanges(seqnames=reads[[1]][vld], ranges=IRanges(start=reads[[3]][vld], width=reads[[4]][vld]), strand=reads[[2]][vld]);
  seqlengths(gr)<-ref.chr[seqlevels(gr)];
  
  cat('Read in total of', length(reads[[1]]), 'sequenced tags.\n');
  cat(length(gr), 'of them have valid alignment to genome.\n');
  rm(reads); # clean up memory
  
  # get depth per base informatioin
  cov<-coverage(gr);
  cov<-cov[names(ref.chr)[names(ref.chr) %in% names(cov)]];
  #---------------------------------------------------------------------------------------------------------------------------#
  #} #@@@@
  
  #@@@@@@
  #if (!exists('bg.stats')) {
  #--------------------------------------------------------------- Background statistics ------------------------------------------------------------------#
  bg.stats<-GetBgStats(bam.file, map.bam.chr.rev[as.vector(seqnames(one.percent))], start(one.percent), end(one.percent), as.vector(strand(one.percent)));
  #--------------------------------------------------------------------------------------------------------------------------------------------------------#  
  if (verbose) cat("Generated genome statistics.\n");
  
  # summarize base frequency
  ex<-colSums(genome$one.percent$base.frequency[bg.stats[[1]], ]); # expected base frequency
  ex<-100*ex/sum(ex);
  ob<-colSums(bg.stats$base$frequency)[names(ex)]; # observed based frequency
  ob<-100*ob/sum(ob);   
  bg.stats$bases$total<-rbind("Expected (%)"=round(ex, 2), "Observed (%)"=round(ob, 2), "Observed/Expected"=round(ob/ex, 2));
  
  # summarize frequency of the first two base combinations
  di<-bg.stats$bases$di.base;
  di.ex<-colSums(genome$one.percent$di.base.frequency[bg.stats[[1]], ]);  # expected frequency of di-base combinations
  di.ex<-di.ex[names(di)];
  di<-t(matrix((di/sum(di))/(di.ex/sum(di.ex)), nc=4)); # normalized to expected frequency
  rownames(di)<-colnames(di)<-c('A', 'C', 'G', 'T');
  names(dimnames(di))<-c('Fisrt base', 'Second base');
  bg.stats$bases$di.normalized<-di;
  
  
  #--------------------------------------------------------------- Whole genome statistics ------------------------------------------------------------------#
  # General
  genome.stats<-list();
  genome.stats$global<-list(
    "Sample name" = sample.name,
    "Genome name" = genome.name,
    "Genome size (bp)" = genome.size,
    "Effective size (bp)" = genome.size - gap.size, 
    "Read length (bp)" = max(width(gr)),
    "Aligned reads" = length(gr),
    "Average depth" = round(weighted.mean(sapply(cov, function(c) mean(as.numeric(c))), as.numeric(sapply(cov, length))), 2),
    "Average sequencing quality" = round(weighted.mean(bg.stats$quality$summary[,2], bg.stats$quality$summary[,4], na.rm=TRUE), 2),
    "Average mapping quality" = round(weighted.mean(as.numeric(names(bg.stats$mapq)), bg.stats$mapq, na.rm=TRUE),  2)
    );
  
  # Percentage of bases over given cutoffs
  depth.counts<-sapply(cov, function(cov) apply(depth.cutoffs, 1, function(c) sum(as.numeric(runLength(cov)[runValue(cov)>=c[1]&runValue(cov)<=c[2]]))));
  depth.counts<-rowSums(depth.counts);
  depth.counts[1]<-max(0, depth.counts[1]-gap.size); # exclude assembly gaps
  depth.percents<-round(depth.counts/(genome.size-gap.size)*100, 2);
  names(depth.counts)<-paste(rownames(depth.cutoffs), ', (', depth.percents, '%)', sep='');  
  genome.stats$depth.counts<-depth.counts;
  
  genome.stats$one.percent<-bg.stats;
  #} #@@@@@@
  #----------------------------------------------------------------------------------------------------------------------------------------------------------#

  t<-targets;
  targets<-targets[as.vector(seqnames(targets)) %in% names(cov)]; # exclude targets whose chromosome names are not found in reference genome
  targets<-targets[as.vector(seqnames(targets)) %in% names(map.bam.chr.rev[!is.na(map.bam.chr.rev)])]; # exclude targets whose chromosome names are not found in bam file
  targets.out<-t[!(names(t) %in% names(targets))]; # Excluded targets
  
  toi<-names(targets)[elementMetadata(targets)$toi==1];
  if (length(toi)==0) toi<-names(targets);
  
  #@@@@@@
  #if (!exists('target.stats')) {  
  #------------------------------------------------------ Statistics of targeted regions ---------------------------------------------------------#
  target.stats<-GetRegionsStats(bam.file, map.bam.chr.rev[as.vector(seqnames(targets))], start(targets), end(targets), as.vector(strand(targets)), names(targets));
  #-----------------------------------------------------------------------------------------------------------------------------------------------#
  if (verbose) cat("Generated target statistics.\n");
  
  #------------------------------------------ Target coverage --------------------------------------------
  # get coverage at each targeted regions and split the targets
  target.stats$coverage$all<-lapply(1:length(targets), function(i) (cov[[as.vector(seqnames(targets))[i]]][start(targets)[i]:end(targets)[i]]));
  names(target.stats$coverage$all)<-names(targets);
  if (verbose) cat('Obtained coverage at targets.\n');

  target.stats$coverage$length<-sapply(target.stats$coverage$all, length);
  target.stats$coverage$mean<-sapply(target.stats$coverage$all, mean);
  target.stats$coverage$cutoff.counts<-sapply(target.stats$coverage$all, function(cov) 
    apply(depth.cutoffs, 1, function(c) sum(as.numeric(runLength(cov)[runValue(cov)>=c[1]&runValue(cov)<=c[2]]))));
  
  depth.counts<-rowSums(target.stats$coverage$cutoff.counts);
  depth.percents<-round(depth.counts/sum(depth.counts)*100, 2);
  names(depth.counts)<-paste(rownames(depth.cutoffs), ', (', depth.percents, '%)', sep='');
  target.stats$coverage$depth.counts<-depth.counts[depth.counts>0];
  
  target.stats$coverage$accu<-quantile(unlist(lapply(target.stats$coverage$all, as.vector), use.name=FALSE), probs=seq(0,1, 0.01));
  
  if (length(targets)>=25000) ce.ind<-sample(1:length(targets), 25000) else ce.ind<-1:length(targets); # to reduce time
  cen.end<-sapply(target.stats$coverage$all[ce.ind], function(x) {
    ws<-ceiling(length(x)/100);                                                             
    stt<-seq(1, length(x), length.out=101);
    end<-pmin(length(x), stt+ws);
    aggregate(x, start=stt, end=end, FUN=mean);
    });

  cen.end<-apply(cen.end, 1, function(a) a/target.stats$coverage$mean[ce.ind]);
  target.stats$coverage$cen.end<-cen.end[!is.na(target.stats$coverage$mean[ce.ind]), ];
  if (verbose) cat('Obtained center vs. end relative depth.\n');  

  target.stats$stats$depth<-sapply(target.stats$coverage$all, function(x) summary(as.numeric(x))[1:6]);
  target.stats$stats$qual<-sapply(target.stats$qual, function(x) summary(as.numeric(x))[1:6]);
  target.stats$stats$mapq<-sapply(target.stats$mapq, function(x) summary(as.numeric(x))[1:6]);
  
  #------------------------------------------------------------------------------------------------------
  
  target.stats$basic<-list(
    "Number of targets" = length(targets),
    "Missing targets" = length(targets.out),
    "Targets of interest" = length(toi),
    "Total size of targets (bp)" = sum(as.numeric(width(targets))),
    "Total on-target reads" = sum(as.numeric(target.stats$count)),
    "Average depth" = round(weighted.mean(sapply(target.stats$coverage$all, function(c) mean(as.numeric(c))), as.numeric(sapply(target.stats$coverage$all, length))), 2),
    "Average sequencing quality" = round(mean(unlist(target.stats$qual, use.names=FALSE), na.rm=T), 2),
    "Average mapping quality" = round(mean(unlist(target.stats$mapq, use.names=FALSE), na.rm=T), 2)    
    );
  #}#@@@@@

  #------------------------------- Run Sweave ----------------------------------------------------------#
  Sweave('bamchop.rnw');
  
  system('pdflatex bamchop.tex');
  system('pdflatex bamchop.tex');
  #-----------------------------------------------------------------------------------------------------#
  
  # move files
  folder.name<-paste(path.out, sample.name, sep='/');
  if (!file.exists(folder.name)) dir.create(folder.name);
  file.copy('bamchop.pdf', paste(path.out, '/', sample.name, '.pdf', sep=''));
  pdfs<-dir();
  pdfs<-pdfs[grep('.pdf$', pdfs)];
  pdfs<-pdfs[grep('^bamchop', pdfs)];
  file.rename(pdfs, paste(folder.name, pdfs, sep='/'));
  file.remove(paste('bamchop.', c('aux', 'log', 'out', 'tex', 'toc'), sep=''));  
}

########### TODO: NEED REVISION
#------------------------------- WRITE OUT ----------------------------------------------------------#
if (write.out) {  
  # UCSC file, depth at targeted regions, background is set to 0
  c<-coverage(targets);
  for (i in 1:length(c)) if (max(c[[i]])>1) c[[i]][c[[i]]>1]<-1; # In case there are overlapping targets
  Depth2BedGraph(c*cov[names(c)], sample.name, folder.name, 'targets', db='hg19')
  
  # position with higher depth than cutoff within the background
  cov.bg<-cov; 
  cov.bg[names(c)]<-(1-c)*cov.bg[names(c)];
  bg.cut<-quantile(target.stats$coverage$mean)[2];
  cov.bg[cov.bg<bg.cut]<-0;
  Depth2BedGraph(cov.bg, sample.name, folder.name, 'high_bg', db='hg19');

  save(gr, file=paste(folder.name, '/gr.rdata', sep=''));
  save(cov, file=paste(folder.name, '/cov.rdata', sep=''));  
}
#----------------------------------------------------------------------------------------------------#

dates[[length(dates)+1]]<-date();

#Backup the parameter file
file.copy('bamchop_param.r', paste(path.out, '/bamchop_param.r', sep=''));