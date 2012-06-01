# Given a list of Rle objects with each element corresponding to a chromosome and it read depth at each base, draw a landscape plot to visualize the depth information altogether
# Optional, highlight a given list of regions of interest
PlotLandscape<-function(cov, roi=NA, ws=10, step=10, same.x=TRUE, same.y=TRUE, 
                        y.scale=list('3SD', 'MAX', '99.9%', 'LOG', 'CONSTANT'=600), 
                        col.bg='#DDDDDD', col.fg='steelblue', col.roi='orange', col.chr='white',
                        pdf=FALSE) {

 
  len<-sapply(cov, length); # length of chromosomes, in million bp
  ws<-ws*1000; # window size to take averages
  ws<-max(max(len)/100000, ws); # minimal widnows size to reduce calculation time
  step<-step*1000; # step size to take averages
  step<-max(max(len)/100000, step); # minimal step size to reduce plotting time

  # calculate XY coordinates, X is the chromosome location (in million bp) and Y is mean depth on a chromosome
  getXY<-function(c, ws, step) { 
    start<-seq(1, length(c), step);
    end<-pmin(length(c), start+ws-1); 
    list(X=start/2000000+end/2000000, Y=aggregate(c, start=start, end=end, FUN=mean));
  }
  
  XY<-lapply(cov, function(cov) getXY(cov, ws, step));
  # x-axis limit of each sub-plot
  xlim<-len/10^6;
  if (same.x)   xlim<-rep(max(xlim), length(cov));
  
  # y-axis limit of each sub-plot
  {if ((class(y.scale)=='numeric'|class(y.scale)=='integer') & y.scale[[1]]>0) ylim<-rep(y.scale, length(XY)) # specified number 
  else if (identical(y.scale[[1]], 'MAX')) ylim<-sapply(XY, function(x) max(x[[2]])) # the highest
  else if (identical(y.scale[[1]], '99.9%')) ylim<-sapply(XY, function(x) quantile(x[[2]], probs=seq(0, 1, 0.001)[1000])) # 99.9 percentile
  else if (identical(y.scale[[1]], 'LOG')) ylim<-log10(sapply(XY, function(x) max(x[[2]]))) # log10 of the highest
  else ylim<-sapply(XY, function(x) mean(x[[2]])+3*sd(x[[2]]));
   }
  if (same.y) ylim<-rep(max(ylim, na.rm=TRUE), length(XY));

  # chromosome names to be labeled on the plot
  chr.names<-names(cov);
  chr.names<-sapply(1:length(cov), function(i) if (is.null(chr.names[i])) i 
                    else if (is.na(chr.names[i])) i 
                    else chr.names[i]);
  chr.names<-sub('chr', '', chr.names, ignore.case=TRUE);
  
  # plot settings
  if (pdf) pdf('landscape.pdf', w=8.5, h=11); # A4 size
  par(mar=c(0,3,0,0.5), oma=c(0.5,0,0.5,0), mfrow=c(length(cov), 1), bg=col.bg);
  
  # text size of chromosome names
  chr.h<-strheight(chr.names, units='figure', cex=1);
  chr.w<-20*strwidth(chr.names, units='figure', cex=1);
  
  # return values
  out<-list(all=list(X=list(), Y=list()));
  if (class(roi)=='GRanges') out$roi<-list(X=list(), Y=list());
  
  # plot chromosomes one by one
  for (i in 1:length(XY)) {
    plot(0, type='n', xlim=c(0, 1.05*xlim[i]), ylim=c(0, ylim[i]), axes=FALSE, xaxs='i', yaxs='i');
    
    # grids
    grids<-seq(0, xlim[i], 10); 
    abline(v=grids, lty=1, lwd=1, col='white'); 
    abline(v=xlim[i], lty=1, lwd=1.5, col='black');
  
    # Chromosome names
    text(1.025*xlim[i], ylim[i]/2, label=chr.names[i], col=col.chr, cex=0.75/max(c(chr.h, chr.w))); 
    
    x<-XY[[i]][[1]]; 
    y<-XY[[i]][[2]];
    if(identical(y.scale, 'LOG')) y<-log10(pmax(1, y)); # convert to log scale
  
    # plot y-axis
    if (identical(y.scale, 'LOG')) {
      ytick<-0:floor(ylim[i]);
      ylab<-10^ytick;
    }
    else {
      ytick<-pretty(c(0, ylim[i]), n=5);
      ylab<-ytick<-ytick[-length(ytick)];
    }
    axis(2, las=2, at=ytick, labels=ylab, cex.axis=2/3);
                  
    polygon(c(0, x, 0), c(0, y, 0), col=col.fg, border=col.fg);
    
    out[[1]]$X[[names(cov)[i]]]<-x;
    out[[1]]$Y[[names(cov)[i]]]<-y;
    
    # plot regions of interest if there is any on the chromosome
    if (class(roi)=='GRanges') {  
      gr<-roi[seqnames(roi)==names(cov)[i]]; 
      if (length(gr)>0) {
       stts<-floor(start(gr)/ws);
       ends<-ceiling(end(gr)/ws);
       ind<-unlist(sapply(1:length(stts), function(i) stts[i]:ends[i]));
       y[-ind]<-0;
       polygon(c(0, x, 0), c(0, y, 0), col=col.roi, border=col.roi);
       
       out[[2]]$X[[names(cov)[i]]]<-x;
       out[[2]]$Y[[names(cov)[i]]]<-y;
      }
      else {
        out[[2]]$X[[names(cov)[i]]]<-x;
        out[[2]]$Y[[names(cov)[i]]]<-rep(0, length(y));
      }
    }
    
    box(lwd=1.5, col='black');
  }
  
  if (pdf) dev.off();
  
  out
} 
