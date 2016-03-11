\name{production_url_list}
\alias{production_url_list}
\docType{data}
\title{URLs of NASA GPM IMERG data}
\description{This is an example of the NASA GPM IMERG v3 production URLs as obtained using the rcurl_download function of this package. The files contain 49hrs of global files from 2014-4-01 to 2014-04-02.}
\usage{data("production_url_list")}
\format{
  A data frame with 49 observations on the following variable.
  \describe{
    \item{\code{V1}}{a factor with levels defining the URLs from the NASA GPM IMERG ftp server for IMERG production data}
  }
}
\details{This is an example of the NASA GPM IMERG v3 production URLs as obtained using the rcurl_download function of this package. The files contain 49hrs of global files from 2014-4-01 to 2014-04-02.}
\source{NASA PPS, http://pmm.nasa.gov/data-access/downloads/gpm}
\examples{
data(production_url_list)
str(production_url_list)
production_url_list[1,1]
}
\keyword{datasets}
