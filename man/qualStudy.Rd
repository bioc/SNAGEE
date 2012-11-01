\name{qualStudy}
\alias{qualStudy}
\title{Quality of a study}
\description{
	Calculate the quality of a study.
}
\usage{
qualStudy(d,mode="complete",cc=NULL,disattenuate=TRUE)
}
\arguments{
	\item{d}{The study data. If an Eset, \link{toSnageeFormat} is called on it.
	    Otherwise, must be a list with fields 'genes' containing the
		vector of gene IDs (from NCBI's Gene DB) and 'data' containing the actual data. }
	\item{mode}{Which gene-gene correlation matrix should be used.
		Can be 'complete' (using all platforms) or 'woAffy' (without the Affy platforms). }
	\item{cc}{Can be used if wishing to use a custom gene-gene correlation matrix.
		Must be a list with fields 'g' containing the gene IDs and 'cc' containing the
		  (upper triangular part of the) correlations. }
	\item{disattenuate}{Should the qualities be disattenuated?}

}
\details{
	The function calculates the quality of a study. The numerical values of the study
	(the 'data' field) should be in log-scale, and normalized. It is recommended to used medpolish
	on the data.
	
	Each gene should only appear once in the gene list. Duplicated genes must be merged before
	using the function.
	
	The \code{mode} \code{'woAffy'} may be useful to compare Affymetrix to not Affymetrix studies.
	As the median gene correlation matrix was calculated with a majority of Affymetrix platforms, those
	platforms tend to be given higher quality than the others with the \code{'complete'} mode,
	which may be misleading. 
}
\seealso{\link{SNAGEE}, \link{qualSample}, link{toSnageeFormat}
}
\examples{
# Get the list of genes
geneList = getCC()$g;
# Create a random data set
d=list(genes=geneList, data=matrix(rnorm(length(geneList)*50),ncol=50));
# And calculate its quality (it's going to be close to 0)
qualStudy(d, disattenuate=FALSE);
}
