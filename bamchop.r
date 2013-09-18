########################### Mandatory variables ##########################################################################
bamchop.path<-"/home/zhangz/hts/scripts/bamchop.1.01";

genome.names<-c('ce10');

bamfiles<-c(
'ACAGTG_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/ACAGTG_1Aligned.out.sorted.bam", 
'ACAGTG_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/ACAGTG_2Aligned.out.sorted.bam", 
'CAGATC_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CAGATC_1Aligned.out.sorted.bam", 
'CAGATC_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CAGATC_2Aligned.out.sorted.bam", 
'CCGTCC_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CCGTCC_1Aligned.out.sorted.bam", 
'CCGTCC_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CCGTCC_2Aligned.out.sorted.bam", 
'CGATGT_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CGATGT_1Aligned.out.sorted.bam", 
'CGATGT_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CGATGT_2Aligned.out.sorted.bam", 
'CTTGTA_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CTTGTA_1Aligned.out.sorted.bam", 
'CTTGTA_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/CTTGTA_2Aligned.out.sorted.bam", 
'GCCAAT_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/GCCAAT_1Aligned.out.sorted.bam", 
'GCCAAT_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/GCCAAT_2Aligned.out.sorted.bam", 
'GTGAAA_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/GTGAAA_1Aligned.out.sorted.bam", 
'GTGAAA_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/GTGAAA_2Aligned.out.sorted.bam", 
'TGACCA_1'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/TGACCA_1Aligned.out.sorted.bam", 
'TGACCA_2'="/home/zhangz/hts/projects/mf001/2013-06_RNAseq_Pilot/10k_run_single_ended/alignments/star/TGACCA_2Aligned.out.sorted.bam"
  
);
################################################################################################################

##################### EXTRA SECTIONS IN REPORT ##########################
# Application type of the sequencing data, such as RNA-seq and Exome-seq
# A corresponding section will be added into the report if specified
type<-'';
#########################################################################
roi<-NA; # Regions of interest; GRanges object or full name of a BED file
soi<-NA; # SNPs of interest; GRanges object or full name of a BED file
#########################################################################

############ Extra parameters ############ 
extra.param<-list(bamchop.path=bamchop.path);

extra.param$project.name<-'Worm RNA-seq 10K one end';
extra.param$prepared.by<-'Zhe Zhang';
extra.param$affiliation<-'Center for Biomedical Informatics at CHOP';

# Location of image with the icon to be put in footer
extra.param$logo<-paste(bamchop.path, 'logo.jpg', sep='/');
################################################################################################################

# Landscape plotting 
extra.param$landscape.ws<-100;  # Window size of the regions to calculate average depth, in unit of 1kb, default is 10kb
extra.param$landscape.type<-.5; 

# estimated fragment length
extra.param$fragment.length<-200;

# Minimum chromosome size
extra.param$include.chr<-10^7;

################################################################################################################
################################################################################################################
debug<-FALSE; # if true, make a fast test run without re-generating all statistics
rerun<-FALSE; # if true, load pre-existing bamchop object to generate a report

# loading custome functions
source.path<-paste(bamchop.path, '/source', sep='');
sources<-dir(source.path);
sources<-sources[grep('R$', sources)];
sapply(paste(source.path, sources, sep='/'), source)->x;

library(xtable);

# prepare Sweave template
file.copy(extra.param$logo, 'logo.jpg');
file.copy(paste(bamchop.path, 'Sweave.sty', sep='/'), 'Sweave.sty');
CreateRnw(bamchop.path, type, roi, soi);

# make sure each bamfile is associated with a corresponding genome name
if (length(genome.names)!=length(bamfiles)) genome.names<-rep(genome.names[1], length(bamfiles));

# make sure each bamfile has a unique name
nm<-as.vector(sub('.bam', '', sapply(strsplit(bamfiles, '/'), function(x) (x[length(x)]))));
while(length(nm[duplicated(nm)])>0) {
  nm[duplicated(nm)]<-paste(nm[duplicated(nm)], 'dup', sep='_');
}
if (is.null(names(bamfiles))) names(bamfiles)<-nm  else names(bamfiles)[is.na(names(bamfiles))|names(bamfiles)=='']<-nm[is.na(names(bamfiles))|names(bamfiles)==''];

########################################################################################################################

err<-list();
for (i in 1:length(bamfiles)) { print(bamfiles[i]);
  err[[names(bamfiles)[i]]]<-tryCatch({
      if (rerun) load(paste(names(bamfiles)[i], '/', names(bamfiles)[i], '.rdata', sep=''));
    
      if (!rerun&!debug) {
        genome<-LoadGenome(genome.names[i], bamchop.path);
        if (identical(genome, NA)) stop(paste('Genome name:', genome.names[i], 'not recognized.'));
      
        bamchop<-SummarizeGenome(bamfiles[i], genome, include.chr=extra.param$include.chr); 
        bamchop$extra.param<-extra.param;
      
        if (type=='chip') bamchop<-SummarizeChip(bamchop);
        
        bamchop$alerts<-CreateAlerts(bamchop);
      }
      
      
      Sweave('bamchop.rnw');
      system('pdflatex bamchop.tex'); 
      system('pdflatex bamchop.tex');
      
      # write outputs
      if (!debug) {
        if (!rerun) if (!file.exists(names(bamfiles)[i])) dir.create(names(bamfiles)[i]);
        
        file.rename('bamchop.pdf', paste(names(bamfiles)[i], '.pdf', sep=''));
        if (!rerun) save(bamchop, file=paste(names(bamfiles)[i], '/', names(bamfiles)[i], '.rdata', sep=''));
      
        bamchop$mapping.location<-NA;
        bamchop$coverage<-NA;
        bamchop$random.subset<-NA;
        if (!rerun) save(bamchop, file=paste(names(bamfiles)[i], '/', names(bamfiles)[i], '_slim.rdata', sep=''));
      
        f<-dir();
        f<-f[grep('.pdf$', f)];
        file.rename('bamchop.tex', paste(names(bamfiles)[i], '/bamchop.tex', sep=''));
        file.rename(f, paste(names(bamfiles)[i], '/', f, sep=''));
      }
  }, error = function(e) e);
}

print(err); 
save(err, file='errors.rdata');

