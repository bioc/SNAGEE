\name{toSnageeFormat}
\alias{toSnageeFormat}
\title{Turns an Eset into a list}
\description{
	Turns an Eset into a list usable by SNAGEE.
}
\usage{
toSnageeFormat(data)
}
\arguments{
	\item{data}{An Eset. If already a list, leaves it as it is.}
}
\details{
	The function turns an Eset into a list usable by SNAGEE. Gene ID annotations are found
	using the annotation slot of the Eset, and the related annotation DB.
	If no annotation DB can be found, gives an error.
	
	In addition, features with identical gene IDs are averaged, and the data are medpolished.
	}
\seealso{\link{SNAGEE}, \link{qualStudy}, \link{qualSample}
}
\examples{
# Get the list of genes
geneList = getCC()$g;
# Create a random data set
d=list(genes=geneList, data=matrix(rnorm(length(geneList)*50),ncol=50));
# And calculate its quality (it's going to be close to 0)
qualStudy(d, disattenuate=FALSE);
}
