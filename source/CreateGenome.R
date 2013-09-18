CreateGenome<-function(BSgenome) {
	library(GenomicFeatures);
	
	genome<-list(name=BSgenome@provider_version, species=BSgenome@species, provider=BSgenome@provider)
	
	seqinfo<-BSgenome@seqinfo;
	seqnames<-seqinfo@seqnames;
	
	# sizes
	len<-seqinfo@seqlengths;	
	names(len)<-seqnames;
	masks<-lapply(seqnames, function(s) {x<-masks(BSgenome[[s]]); if (is.null(x)) Mask(0, start=NULL, end=NULL, width=NULL) else x;});
	w<-sapply(masks, function(x) maskedwidth(x)[1]);
	base<-sapply(seqnames, function(s) alphabetFrequency(BSgenome[[s]]));
	base<-t(base[1:4,]);
	genome$size<-data.frame(Total=len, Effective=len-w, base);
	
	#gaps
	gaps<-lapply(masks, as.data.frame);
	if (length(gaps)==0) genome$gaps<-GRanges() else {
    gaps<-lapply(gaps, function(x) x[x[[1]]=='AGAPS', 2:3]);
	  chr<-rep(seqnames,sapply(gaps, nrow));
	  gaps<-do.call('rbind', gaps);
	  genome$gaps<-GRanges(chr, IRanges(gaps[,1], gaps[,2]), seqlengths=len);
	  names(genome$gaps)<-1:length(genome$gaps);
	}
  
	# Get refseq from UCSC
	if (genome$provider=='UCSC') {
		db<-makeTranscriptDbFromUCSC(genome=genome$name, tablename='refGene');
		saveFeatures(db, paste(genome$name, '.sqlite', sep=''));
		
		tx<-transcripts(db);
		names(tx)<-elementMetadata(tx)[[1]];
		
		#exon<-exons(db, columns=c("exon_id", 'tx_id'));
		#names(exon)<-elementMetadata(exon)[[1]];
		
		tss<-resize(tx, 1);
		
		promoter<-shift(resize(tss, 1000, fix='end'), 1);
		
		exon<-unlist(exonsBy(db, by=c("tx")));
		cds<-unlist(cdsBy(db, by=c("tx")));
		intron<-unlist(intronsByTranscript(db));
		u5<-unlist(fiveUTRsByTranscript(db));
		u3<-unlist(threeUTRsByTranscript(db));
		
		genome$tss<-tss;
		genome$features<-list(Promoter=promoter, '5-UTR'=u5, Coding=cds, Intron=intron, '3-UTR'=u3, Exon=exon);
		
		cov<-lapply(genome$features, coverage);
    genome$total.feature.length<-sapply(cov, function(c) sum(sapply(c, function(c) length(c[c>0]))));
	}
	
	genome;
}