# Parse cigar string
ParseCigar<-function(cigar, code=NA) {
# cigar   cigar string
# code    if NA, parse all code; otherwise, parse given code
  
  x<-as.integer(strsplit(cigar, '[a-z|A-Z]')[[1]]);
  y<-strsplit(cigar, '[0-9]')[[1]];
  names(x)<-y[y!=''];
  
  if (identical(NA, code)) x else x[code];
}