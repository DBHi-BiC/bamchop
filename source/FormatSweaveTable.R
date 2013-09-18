# Format a table in Sweave document with default parameters
FormatSweaveTable<-function(t, row=NA, col=TRUE, caption=NULL, longtable=FALSE, digits=2, scientific=FALSE, force.size=NA) {
# row   if equals to "", the row names will not be included; otherwise, the row names will be included at the first column with the given value as column name
	library(xtable);
  
  t<-as.data.frame(t);
	
  # if row name should be included, set it as the first column
  if (!is.na(row[1])) {
    t<-cbind(rownames(t), t);
    names(t)[1]<-row;
  }
  
  # format columns
	for (i in 1:ncol(t)) {		
		if (class(t[[i]])=='integer') t[[i]]<-format(t[[i]], big.mark=',')
		else if (class(t[[i]])=='numeric') {
      x<-t[[i]];
      if (max(x)==0 & min(x)==0) t[[i]]<-rep('0.00', length(x))
      else if (max(abs(x))/min(abs(x))>=1000 & max(abs(x))>1000) t[[i]]<-sapply(x, function(x) format(x, big.mark=',', digit=max(2, ceiling(log10(x))+2), scientific=FALSE))
      else t[[i]]<-format(x, big.mark=',', digits=max(2, ceiling(log10(min(abs(x[x!=0]))))+2), scientific=scientific);
		}
	}
	
  # vertical line
  align<-paste(c(rep('|r', ncol(t)+1), '|'), collapse='');

  # longtable?
  if (longtable) {
    l<-'longtable';
    size<-'\\tiny';
    ind1<-1;
  }
  else {
    l<-'tabular';
    size<-NULL;
    ind1<-0;
  }
  
  if (!identical(force.size, NA) & class(force.size[1])=='character') size<-force.size[1];

	print(xtable(t, align=align, digits=digits, caption=caption), size=size, floating=FALSE, include.rownames=FALSE, include.colnames=col, tabular.environment=l, caption.placement="top", 
        add.to.row=list(as.list(seq(ind1, nrow(t)-1, 2)), rep("\\rowcolor[gray]{0.9}", length(seq(ind1, nrow(t)-1, 2)))));
	
}
