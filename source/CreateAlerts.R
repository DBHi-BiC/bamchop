CreateAlerts<-function(bamchop) {
  # cutoff values
  selected<-0.01;
  unmapped<-0.5;
  qual<-25;
  mapq<-30;
  dup<-0.8;
  n.base<-0.005;
  gc<-1.1;
  kmer<-0;
  
  stats<-bamchop$stats; 
  summ<-stats$summary;
  
  alerts<-c();

  # percentage of unmapped reads
  pct<-summ['Total upmapped reads']/(summ['Total mapped reads']+summ['Total upmapped reads']);
  if (pct>unmapped) {
    pct<-paste(format(100*pct, digits=2), '%', sep='');
    pct0<-paste(format(100*unmapped, digits=2), '%', sep='');
    alerts<-c(alerts, paste('More than ', pct0, ' (', pct, ') of the total reads in the BAM file were not mapped', sep=''));
  }
  
  # percentage of randomly selected reads
  pct<-length(bamchop$random.subset[['qname']])/summ['Total mapped reads'];
  if (pct<selected) {
    pct<-paste(format(100*pct, digits=2), '%', sep='');
    pct0<-paste(format(100*selected, digits=2), '%', sep='');
    alerts<-c(alerts, paste('Less than ', pct0, ' (', pct, ') of the total reads were randomly selected to summarize', 
                            ' sequencing quality, mapping quality, mismatch frequency and base frequency.', sep=''))
  }
  
  # average sequencing quality
  q<-summ['Average sequencing quality'];
  if (q<qual) alerts<-c(alerts, paste('The average sequencing quality is less than ', qual, ' (', q, ').', sep=''));
  
  # average mapping quality
  q<-summ['Average mapping quality'];
  if (q<mapq) alerts<-c(alerts, paste('The average mapping quality is less than ', mapq, ' (', q, ').', sep=''));
  
  # duplicated mapping
  pct<-1-stats$mapping$duplicates[1,'Percentage']/100;
  if (pct>dup) {
    pct<-paste(format(100*pct, digits=2), '%', sep='');
    pct0<-paste(format(100*dup, digits=2), '%', sep='');
    alerts<-c(alerts, paste('More than ', pct0, ' (', pct, ') of the mapped reads share the same mapping locations with other reads (duplication).', sep=''));    
  }
  
  # GC contents
  rt<-bamchop$stats$bases$observed.vs.expected["Observed/Expected(%)", 'GC']/100;
  if (rt>gc) {
    pct<-paste(format(100*rt, digits=4), '%', sep='');
    pct0<-paste(format(100*gc, digits=4), '%', sep='');
    alerts<-c(alerts, paste('GC contents of reads are more than ', pct0, ' (', pct, ') of expected percentage.', sep=''));    
  }
  
  # percentage of N base
  pct<-summ['Base N%']/100;
  if (pct>n.base) {
    pct<-paste(format(100*pct, digits=2), '%', sep='');
    pct0<-paste(format(100*n.base, digits=2), '%', sep='');
    alerts<-c(alerts, paste('More than ', pct0, ' (', pct, ') of the bases are Ns.', sep=''));    
  }
  
  # kmer outliers
  k<-stats$bases$position.specific$kmer$k;
  o<-stats$bases$position.specific$kmer$outliers;
  if (nrow(o)>0) alerts<-c(alerts, paste('There are ', nrow(o), ' ', k, '-mers overrepresented at either end of reads.', sep=''));
  
  alerts;
}