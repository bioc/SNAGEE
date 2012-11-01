qualStudy = function(d,mode="complete",cc=NULL,disattenuate=TRUE)
{  d=toSnageeFormat(d);

   if (is.null(cc)) { cc = getCC(mode=mode); }
	
   w = (d$genes %in% cc$g) & apply(is.finite(d$data),1,sum)>.75*ncol(d$data);
   if (sum(w)<20) { stop("Not enough probes\n"); }	
   d$data=d$data[w,]; d$genes=d$genes[w];
  
   if (any(table(d$genes)>1)) { stop("Each gene should only appear once"); } 
   w = match(cc$g,d$genes);
   cb = cor(t(d$data[w,]),use="pairwise.complete.obs");
	
   if (disattenuate)
   {	cb = cb[upper.tri(cb)]
   	    cb = atanh(cb);
		cb[!is.finite(cb)]=NA;
		N = ncol(d$data); er = 1/(N-1)+2/((N-1)^2);
		s = sd(cb,na.rm=TRUE);
		return(cor(cb,atanh(cc$cc),use="p")*sqrt((s^2)/(s^2-er)));
	}
	else 
	{ return(cor(cb[upper.tri(cb)],cc$cc,use="p")); }
}

qualSample = function(data, mode="complete", cc=NULL, multicore=FALSE)
{	if (multicore && !require("multicore"))
	{ warning("multicore not installed. Using only one core.\n"); multicore=FALSE; }
	
	data=toSnageeFormat(data);
	
	if (is.null(cc)) { cc = getCC(mode=mode); }

	if (any(!is.finite(data$data)))
	{ stop("Does not work with missing data. Use impute or remove them first.") } 
	w = (data$genes %in% cc$g);
	data$data=data$data[w,]; data$genes=data$genes[w];
	if (any(table(data$genes)>1)) { stop("Each gene should only appear once"); }

	#message("Calculations are on ", nrow(data$data), " genes.\n");
	if (nrow(data$data)<20) { stop("++ Not enough genes ++\n\n"); }
	
	w = match(data$genes,cc$g);
	t=matrix(NA,nrow=length(cc$g),ncol=length(cc$g));
	t[upper.tri(t)] = cc$cc;
	t[lower.tri(t)] = t(t) [lower.tri(t)]
	t=t[w,w];
	ut=upper.tri(t);
	t=t[ut];
	
	data=data$data;
	n=ncol(data);
	sBase = 1/sd(cor(t(data))[ut]);
	
	t=t-mean(t); t=t/sd(t);
	
	s=apply(data,1,sum);
	c=tcrossprod(data);
	si=apply(data^2,1,sum);
	c=(n-1)*c - tcrossprod(s);
	
	calcCor = function(i)
	{   x = data[,i];
        tmp = ( (n-1)*(si-x^2) - (s-x)^2) ^ -.5;
		tmp= ( c- tcrossprod(n*x-s,x) + tcrossprod(x,s) ) * tcrossprod(tmp);
		tmp = tmp[ut];
		
		tmp = tmp-mean(tmp); tmp=tmp*sBase;
		return(mean( (tmp-t)^2 ));
	}

	if (multicore)
	{ res=mclapply(1:ncol(data), calcCor); }
	else
	{ res=lapply(1:ncol(data), calcCor); }
	res = unlist(res);
	
	res = (res - median(res))/mad(res); 
	return(res);
}

toSnageeFormat = function(data)
{  if (is(data, "list")) { return(data); }
   if (!is(data, "eSet")) { stop("Data must be a list or an eSet"); }
   
   dbname = paste(annotation(data), ".db", sep="");
   if (!require(dbname, character.only=TRUE)) { stop("Annotation package ", dbname, " cannot be found"); }
   
   d = exprs(data);
   if (!any(d <= 0) && any(d > 1e3)) # Try to guess if not in log
   { d = log(d); } 
   
   gDB = paste(annotation(data),"ENTREZID", sep="");
   g = eval(parse(text=gDB));
   g = as.numeric(selectMethod("as.list", "AnnDbBimap")(g));
   d = d[!is.na(g),]; g=g[!is.na(g)];
   if (any(table(g) > 1))
   {    cat("Some genes appear more than once - collapsing using the mean.\n");
        y <- rowsum(d, g, reorder = FALSE, na.rm = TRUE) # Taken from avereps in limma
        n <- rowsum(1L - is.na(d), g, reorder = FALSE)
        d = y/n;
        g = unique(g);
   } 
   
   d = medpolish(d, trace.iter = FALSE)$residuals;
   return(list(data=d, genes=g));
}
