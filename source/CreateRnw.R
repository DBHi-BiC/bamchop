CreateRnw<-function(bamchop.path, type, roi, soi) {
  readIn<-function(filename) {
    l<-c();
    con<-file(filename, open="r");
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
      l<-c(l, oneLine);
    }
    close(con);
    l
  }
  
  lines<-readIn(paste(bamchop.path, '/templates/', 'header.rnw', sep=''));
  lines<-append(lines, readIn(paste(bamchop.path, '/templates/', 'main.rnw', sep='')));
  
  if (type=='chip') lines<-append(lines, readIn(paste(bamchop.path, '/templates/', 'chip.rnw', sep='')));
  
  lines<-append(lines, readIn(paste(bamchop.path, '/templates/', 'alerts.rnw', sep='')));
  
  lines<-append(lines, '\\end{document}');
  
  con<-file('bamchop.rnw');
  writeLines(lines, con);
  close(con);
}