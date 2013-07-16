\name{SNAGEE-package}
\alias{SNAGEE-package}
\alias{SNAGEE}
\docType{package}
\title{
Signal-to-Noise Applied to Gene Expression Experiments
}
\description{
Signal-to-Noise Applied to Gene Expression Experiments
}
\details{
\tabular{ll}{
Package: \tab SNAGEE\cr
Version: \tab 0.99.0\cr
Date: \tab 2012-01-26\cr
Depends: \tab R (>= 2.6.0)\cr
Imports: \tab SNAGEEdata\cr
Suggests: \tab ALL\cr 
Enhances: \tab parallel\cr
License: \tab Artistic-2.0\cr
URL: \tab http://fleming.ulb.ac.be/SNAGEE\cr
}

Index:
\preformatted{
qualStudy               Quality of a study
qualSample              Quality of samples in a study
toSnageeFormat          Turns an Eset to a list usable by SNAGEE
}


}
\author{
David Venet <davenet@ulb.ac.be>

Maintainer: David Venet <davenet@ulb.ac.be>
}
\examples{
# Get the list of genes
geneList = getCC()$g;
# Create a random data set
d=list(genes=geneList, data=matrix(rnorm(length(geneList)*50),ncol=50));
# Calculate its quality (it's going to be very close to 0)
qualStudy(d, disattenuate=FALSE);
# Calcuate individual sample qualities
qs = qualSample(d);
}
