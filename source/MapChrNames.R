MapChrNames<-function(q, r) {
	# q		query, vector of integers named by names of query chromosomes, values equal to chromosome lengths
	# r		reference, vector of integers named by names of referece chromosomes, values equal to chromosome lengths
	
	names(q)[is.na(names(q))]<-q[is.na(names(q))];
	
	# first attempt, the same name & the same length
	mp<-sapply(names(q), function(nm) {mp<-names(r)[names(r)==nm&r==q[nm]]; if (length(mp)==0) NA else mp[1];})
	
	# second attempt, the same name with different cases
	if (length(mp[is.na(mp)>0])) {
		q0<-sub('chr', '', tolower(names(q)));
		names(q0)<-names(q);
		r0<-sub('chr', '', tolower(names(r)));
		names(r0)<-names(r);
		mp[is.na(mp)]<-sapply(names(q)[is.na(mp)], function(nm) {mp<-names(r)[r0==q0[nm]&r==q[nm]]; if (length(mp)==0) NA else mp[1];})
	}

	# third attempt, just same name
	if (length(mp[is.na(mp)>0])) {
		mp[is.na(mp)]<-sapply(names(q)[is.na(mp)], function(nm) {mp<-names(r)[names(r)==nm]; if (length(mp)==0) NA else mp[1];})
	}
	
	# fourth attempt, just the same name with different cases
	if (length(mp[is.na(mp)>0])) {
		q0<-sub('chr', '', tolower(names(q)));
		names(q0)<-names(q);
		r0<-sub('chr', '', tolower(names(r)));
		names(r0)<-names(r);
		mp[is.na(mp)]<-sapply(names(q)[is.na(mp)], function(nm) {mp<-names(r)[r0==q0[nm]]; if (length(mp)==0) NA else mp[1];})
	}

	# last attempt, just the same length
	if (length(mp[is.na(mp)>0])) {
		mp[is.na(mp)]<-sapply(names(mp)[is.na(mp)], function(nm) {mp<-names(r)[r==q[nm]]; if (length(mp)==0) NA else mp[1];} )
	}
	
	mp;
}