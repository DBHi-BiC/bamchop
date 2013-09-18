LoadGenome<-function(genome.name, bamchop.path) {
  
  g<-tolower(genome.name); 
  
  if (g %in% c('hg19', 'human', 'homo sapiens', 'hs', 'h.s.', 'hsapiens', 'grch37')) {
    eval(parse(text=load(paste(bamchop.path, '/database/hg19.rdata', sep=''))));
  }
  else if (g %in% c('hg18', 'ncbi36')) {
    eval(parse(text=load(paste(bamchop.path, '/database/hg18.rdata', sep=''))));
  }
  else if (g %in% c('mm9', 'mouse', 'mouse musculus', 'mmusculus', 'ncbi37')) {
    eval(parse(text=load(paste(bamchop.path, '/database/mm9.rdata', sep=''))));
  }
  else if (g %in% c('ce10', 'worm', 'c. elegans', 'celegans', 'ws220')) {
    eval(parse(text=load(paste(bamchop.path, '/database/ce10.rdata', sep=''))));
  }
  else if (g %in% c('mm8')) {
    NA
  }
  else NA;
  
}