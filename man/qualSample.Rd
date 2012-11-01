\name{qualSample}
\alias{qualSample}
\title{Quality of samples in a study}
\description{
	Calculate the relative quality of all samples from a study.
}
\usage{
qualSample(data,mode="complete",cc=NULL,multicore=FALSE)
}
\arguments{
	\item{data}{The study data. If an Eset, \link{toSnageeFormat} is called on it.
	    Otherwise, must be a list with fields 'genes' containing the
		vector of gene IDs (from Entrez) and 'data' containing the gene expression data. }
	\item{mode}{Which gene-gene correlation matrix should be used.
		Can be 'complete' (using all platforms) or 'woAffy' (without the Affy platforms). }
	\item{cc}{Can be used if wishing to use a custom gene-gene correlation matrix.
		Must be a list with fields 'g' containing the gene IDs and 'cc' containing the
		 (upper triangular part of the) correlations. }
	\item{multicore}{Should the parallel version be used? This is based on the multicore package,
		if that package cannot be loaded it will fall back on single core, with a warning.}
}
\details{
	The function calculates the quality of all samples in a study. Lower values are of lower quality.
	The numerical values of the study
	(the 'data' field) should be in log-scale, and normalized. It is recommended to used medpolish
	on the data.
	
	Each gene should only appear once in the gene list. Duplicated genes must be merged before
	using the function. Non-finite values should also be removed first (using the impute package for
	instance).

}
\seealso{\link{SNAGEE}, \link{qualStudy}, \link{toSnageeFormat}
}
\examples{
# Get the list of genes
geneList = getCC()$g;
# Create a random data set
d=list(genes=geneList, data=matrix(rnorm(length(geneList)*50),ncol=50));
# And calculate the quality of the samples (they are all about the same)
qualSample(d);
}
