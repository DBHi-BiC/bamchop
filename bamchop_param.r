# Required parameters
bams<-"/home/zhangz/hts/scripts/bamchop/examples/exome/bams.txt"; 
targets<-"/home/zhangz/hts/scripts/bamchop/examples/exome/exons.txt";
path.out="/home/zhangz/hts/scripts/bamchop/examples/exome";

# Optional parameters
project.name='Whole exome resequencing';
prepared.by='Zhe Zhang';
affiliation='Center for Biomedical Informatics, CHOP';
genome.name='hg19';


# Parameters for advanced users, change with caution
verbose=TRUE;

logo='./logo.jpg';

# mapping chromosome names?
map.chr.names<-TRUE;

# write out files (R objects and full tables)
write.out=FALSE;

# Landscape plotting 
landscape.ws<-10;  # Window size of the regions to calculate average depth, in unit of 1kb, default is 10kb
landscape.type<-'LOG';

# depth cutoffs for pie charts and summary tables
depth.cutoffs<-cbind(min=c(0,1,10,30,100), max=c(0,9,29,99,Inf)); # cutoff values for depth catogories
rownames(depth.cutoffs)<-c('0', '1 - 9', '10 - 29', '30 - 99', '100+');
variant.call.cutoff<-30; # required depth to make a variant call
