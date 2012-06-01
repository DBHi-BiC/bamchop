Plot3DPie<-function (v, main='', border='grey', height=0.02, explode=0.03, radius=3/4) {
  ver<-as.numeric(R.Version()$minor);
  
  par(mfrow=c(1,1), mar=c(0,0,0,0), omi=c(0,0,0,0));
  if (ver>=15) angle<-pie3D(v, theta=pi/6, border=border, height=height, explode=explode, radius=radius)
  else angle<-pie3D(v, theta=pi, border=border, height=height, explode=explode, radius=radius);

  srt<-angle*180/pi;
  srt[srt>90 & srt<180]<-srt[srt>90 & srt<180]+180;
  srt[srt>180 & srt<270]<-srt[srt>180 & srt<270]-180;

  if (ver>=15) sapply(1:length(angle), function(i) pie.labels(0, 0, angle[i], names(v)[i], adj=0.5, radius=min(1, 1.2*radius), border=NA, bg=NA, srt=srt[i], cex=1.5, pos=3))->x
  else sapply(1:length(angle), function(i) pie.labels(0, 0, angle[i], names(v)[i], adj=0.5, radius=0.9*radius, border=NA, bg=NA, srt=srt[i], cex=1.5, pos=3))->x;
  
  title(main=main, cex.main=3, line=-2);
}

PlotLogBar<-function(v, main='', col='orange', space=0, las=3, xlab='', ylab='Percentage (%)') {
  if (is.na(xlab) | xlab[1]=='') mar1<-2 else mar1<-4.5;
  if (is.na(ylab) | ylab[1]=='') mar2<-2 else mar2<-5;
  if (is.na(main) | main[1]=='') mar3<-0.5 else mar3<-3;
  par(mar=c(mar1, mar2, mar3, 1)); 
  
  ticks<-c(seq(0.01, 0.1, 0.01), seq(0.2, 1, 0.1), 2:10, seq(20, 100, 10));
  
  barplot(v, log='y', ylim=c(0.9*min(v), min(100, 1.1*max(v))), xlab=xlab, ylab=ylab, yaxt='n', space=space, las=las, cex.lab=2, col='#FFFFFFFF'); 
  abline(h=ticks, col='lightgrey', lty=2, lwd=0.5); 
  abline(h=10^(-6:2), col='darkgrey'); 
  barplot(v, log='y', col=col, add=TRUE, space=space, names.arg=NA, yaxt='n'); 
  axis(2, at=ticks, labels=gsub(' ', '', format(ticks, drop0trailing=TRUE, scientific=FALSE))); 
  title(main=main, cex.main=3);
  box(lwd=2);
}

# Plot the base frequncy logo
PlotBaseLogo<-function(freq, main='', ic.scale=FALSE, xaxis=FALSE) {
  plot.new();
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  par(mar=c(0, 4, mar3, 1));
  
  seqLogo(freq, ic.scale=ic.scale, xaxis=xaxis);
  title(main=main, cex.main=3);    
}

# Plot overall distribution of quanlity score
PlotQualityDens<-function(dens, main='') {
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  par(mfrow=c(1,1), mar=c(4, 4, mar3, 1), omi=c(0, 0, 0, 0));
   
  plot(0, type='n', xlab='', ylab='', xlim=c(min(dens$x), max(dens$x)), ylim=c(min(0, min(dens$y)), 1.05*max(dens$y)), xaxs='i', yaxs='i'); 
  polygon(dens, border='gold', col='#888888');
  title(ylab='Density', xlab='Quality score', line=2);
  title(main=main, cex.main=3);
  box(lwd=2);
}

# Plot position-specific percentiles of sequencing quality
PlotQualityPerc<-function(qual, main='') {
# qual  matrix, each row is a percentile, each column is a position
  library(gplots);
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  par(mfrow=c(1,1), mar=c(4, 4, mar3, 1), omi=c(0, 0, 0, 0));
  
  ind<-c(0, 1, 2, 5, 10, 25, 50, 75, 100)+1;
  q<-qual[ind,];
  
  col<-colorpanel(nrow(q), '#FFFF00', '#FF4400');
  
  plot(0, type='n', xlab='', ylab='', xlim=c(1, ncol(q)), ylim=c(min(0, min(q)), max(q)), xaxs='i', yaxs='i'); 
  rect(par()$usr[1], par()$usr[3], par()$usr[2], par()$usr[4], col='#DDDDDD', border=NA);
  title(xlab='Position in reads', ylab='Quality score', line=2, cex.lab=1.5);
  #for (i in nrow(q):1) lines(q[i, ], col='#888888');
  
  for (i in nrow(q):1) polygon(c(1, 1:ncol(q), ncol(q)), c(0, q[i, ], 0), col=col[i], border='#FFFFFF', lwd=.5); 

  # add labels
  X<-sapply(2:7, function(i) {d<-smooth(q[i+1,]-q[i-1,]); d[c(1, length(d))]<-0; as.vector(which(d==max(d))[1])});
  Y<-sapply(2:7, function(i) q[i, X[i-1]]);
  text(X, Y, labels=rownames(q)[2:7], cex=1.5);
  mid<-round(ncol(q)/2);
  min.locs<-which(q[1,]==min(q[1,]));
  max.locs<-which(q[nrow(q), ]==max(q[nrow(q), ]));
  min.loc<-min.locs[which(abs(min.locs-mid)==min(abs(min.locs-mid)))][1];
  max.loc<-max.locs[which(abs(max.locs-mid)==min(abs(max.locs-mid)))][1];
  text(c(min.loc, max.loc), c(min(q), max(q)), labels=c('0%', '100%'), cex=1.5, pos=c(3, 1));
  
  # middle line
  lines(qual[51, ], lwd=3, lty=1, col='#2222FF'); 
  
  title(main=main, cex.main=3);
  box(lwd=2);
}

# Plot accumulative distribution
PlotAccuDepth<-function(d, main='') {
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  par(mfrow=c(1,1), mar=c(4, 5, mar3, 1), omi=c(0, 0, 0, 0));

  x.max<-10^(round(log10(d[length(d)-1]))); # limit of x-axis
  
  plot(d, 100:0, type='l', xlim=c(0, x.max), xlab='', ylab='', xaxs='i', yaxs='i');
  abline(v=seq(0, x.max, x.max/20), h=seq(0, 100, 5), lty=1, col='lightgrey');
  polygon(c(0, d, length(d)), c(0, 100:0, 0), col='#00008888');
  lines(d, 100:0, col='#88000088', lwd=2);
  title(xlab='Depth (X)', ylab='Accumulative Percentage (%)', cex.lab=1.5, main=main, cex.main=3);
  box(lwd=2);
}

# Plot end vs. center depth
PlotCenterEnd<-function(cov, ws=21, main='') {
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  par(mfrow=c(1,1), mar=c(4, 5, mar3, 1), omi=c(0, 0, 0, 0));
    
  y<-runmed(colMeans(cov, na.rm=TRUE), min(ws, ncol(cov)));
  plot(0:(length(y)-1), y, type='n', xaxt='n', ylim=c(0, 1.05*max(y)), xlab='', ylab='', xaxs='i', yaxs='i');
  title(xlab='Relative position in targeted regions', ylab='Relative depth', cex.lab=1.5, main=main, cex.main=3);
  abline(v=seq(0, (length(y)-1), round(length(y)-1)/20), col='lightgrey');
  abline(h=1, col='green');
  polygon(c(0, 0:(length(y)-1), length(y)-1), c(0, y, 0), col='#00008888');
  lines(0:(length(y)-1), y, col='red', lwd=2);
  axis(1, at=c(0, (length(y)-1)/2, (length(y)-1)), labels=c('First', 'Center', 'Last'));
  box(lwd=2);
}

# plot density distribution
PlotDensity<-function(d, main='', xlab='', ylab='', col='#90FF90CC', min.x=0) {
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  if (is.na(xlab) | xlab[1]=='') mar1<-1 else mar1<-4;
  if (is.na(ylab) | ylab[1]=='') mar2<-1 else mar2<-4;
  if (is.na(main) | main[1]=='') mar3<-1 else mar3<-4;
  par(mfrow=c(1, 1), mar=c(mar1, mar2, mar3, 1), omi=c(0, 0, 0, 0));
  
  d<-d[d> -Inf & d< Inf & !is.na(d)];
  
  if (length(d)>0) {
  
  dens<-density(d);
  X<-dens$x;
  Y<-dens$y;
  
  plot(0, type='n', xlim=c(max(min.x, min(X)), max(X)), ylim=c(0, 1.05*max(Y)), xlab='', ylab='', xaxs='i', yaxs='i', xaxt='n', yaxt='n');
  x.tick<-axis(1);
  abline(v=x.tick, col='lightgrey', lty=2);
  polygon(c(0, X, 100), c(0, Y, 0), col=col, border='gold', lwd=2);
  title(xlab=xlab, ylab=ylab, cex.lab=1.5, main=main, cex.main=3);
  box(lwd=2);
  } else {
    plot(0, type='0', xlab='', ylab='', axas=FALSE); 
    box();
  }
}


PlotStamps<-function(d, type=c('line', 'pie', 'polygon', 'density'), max.len=1000, col='#FF0000') { 
  if (class(d)=='matrix' | class(d)=='data.frame') l<-nrow(d) else l<-length(d);
  nc<-ceiling(sqrt(l));
  nr<-ceiling(l/nc);
                par(mfrow=c(nr, nc), mar=c(0,0,0,0), omi=c(0,0,0,0));
  
  if (length(type)==0) type<-'line';
  if (is.na(type[1]) | type[1]=='') type<-'line';
  if (!(type[1] %in% c('line', 'pie', 'polygon', 'density'))) type<-'line';
  
  if (is.null(names(d))) names(d)<-paste('toi', 1:length(d), sep='_');

  # use average if too long
  if (class(d)=='list') {
    d<-lapply(d, function(v) {  
      if (length(v) > max.len) {
        ws<-ceiling(length(v)/max.len);
        stt<-seq(1, length(v), length.out=max.len);
        end<-pmin(length(v), stt+ws);
        aggregate(v, start=stt, end=end, FUN=mean);
      } else v;
      })
    max<-max(sapply(d, max));
  } else if (class(d)=='matrix' | class(d)=='data.frame') {
    max<-max(apply(d, 1, max));
  } else {
    max<-max(d);
  }
  
  # plot regions
  for (i in 1:l) {
    if (class(d)=='list') {
      v<-d[[i]];
      nm<-names(d)[i];
    } else if (class(d)=='matrix' | class(d)=='data.frame') {
      v=d[i, ];
      nm<-rownames(d)[i];
    } else {
      v<-d[i]; 
      nm<-names(d)[i];
    }
    v<-v[!is.na(v)];
    v<-as.numeric(v);

    # plotting ----------------------------------------------------------------------
    if (length(v)>0) {
      if (type[1]=='line') {
        plot(0, 0, type='n', ylim=c(0, max), xlim=c(0, length(v)), axes=FALSE);
        lines(v, col=col)
      } else if (type[1]=='polygon') {
        plot(0, 0, type='n', ylim=c(0, 1.05*max), xlim=c(0, length(v)), yaxs='i', axes=FALSE);
        polygon(c(0, 1:length(v), length(v)), c(0, v, 0), col=col, border=NA);        
      } else if (type[1]=='density') {
        dens<-density(v);
        x<-dens$x; 
        y<-dens$y; 
        plot(0, 0, type='n', ylim=c(0, 1.05*max(y)), xlim=c(0, max(x)), yaxs='i', axes=FALSE);
        polygon(c(x[1], x, x[length(x)]), c(0, y, 0), col=col);
      } else if (type[1]=='pie') {
        pie(v, labels=NA, radius=0.9, init.angle=90, clockwise=TRUE, col=col);
      } else plot(0, type='0', xlab='', ylab='', axas=FALSE); 

    } else plot(0, type='0', xlab='', ylab='', axas=FALSE); 
    
    if (!is.na(nm)) title(main=nm, line=-1.25, cex.main=min(3, 10/nc), col.main='#2222FF');    
    box();
    # ------------------------------------------------------------------------------
  }
}

