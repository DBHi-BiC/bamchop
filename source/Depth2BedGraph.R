Depth2BedGraph<-function(cov, sample.name, folder, suffix, db='') {
# cov, coverage data on all chromosomes
# file.name, full path of output file name
# db, genome version
  
  ##############################################################################################
  BedTrackLine<-function(type='', name='', description='', visibility='', color='', itemRgb='', colorByStrand='', 
                         useScore='', group='', priority='', db='', offset='', url='', htmlUrl='') {
    
    values<-c(type, name, description, visibility, color, itemRgb, colorByStrand, useScore, group, priority, db, offset, url, htmlUrl); 
    names<-c('type', 'name', 'description', 'visibility', 'color', 'itemRgb', 'colorByStrand', 'useScore', 'group', 'priority', 'db', offset, 'url', 'htmlUrl');
    names(values)<-names;
    values[!is.na(values)&values!='']->values; 
    values<-sapply(values, function(x) paste("'", x, "'", sep=''));
    
    string<-paste(names(values), '=', values, sep='');
    line<-paste(string, collapse=' ');
    line<-paste('track', line, sep=' ');
    line;
  }
  ##############################################################################################

file.name<-paste(folder, '/', sample.name, '_', suffix, '.bed', sep='');
  
chr<-rep(names(cov), sapply(cov, function(x) length(x@values))); 
start<-unlist(sapply(cov, start)); 
end<-unlist(sapply(cov, end));
value<-unlist(sapply(cov, function(x) x@values))

bed<-data.frame(chr, start, end, value); 


header=BedTrackLine(type='bedGraph', name=sample.name, description=paste('Depth_', sample.name, sep=''), db=db);
write.table(header, file.name, sep='\t', qu=F, row=F, col=F);
write.table(bed,  file.name, sep='\t', qu=F, row=F, col=F, append=T);

system(paste('gzip -f ', file.name)); # zip bedGraph file and delete original

}
########################

