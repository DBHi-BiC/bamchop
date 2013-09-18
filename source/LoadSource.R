# load source files
LoadSource<-function(bamchop.path) {
  source.path<-paste(bamchop.path, '/source', sep='');
  sources<-dir(source.path);
  sources<-sources[grep('R$', sources)];
  sapply(paste(source.path, sources, sep='/'), source)->x;
}