% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CorrelationIndex.R
\name{CorIndex.all}
\alias{CorIndex.all}
\title{Calculated both the non and group-equalised Phi Coefficient of correlation for all habitats.}
\usage{
CorIndex.all(CorIndexVarInput)
}
\arguments{
\item{CorIndexVarInput}{a custom class "CorIndexVar" containing all values from
from only \link[PhiCor]{CorIndex}. It calls \link[PhiCor]{CorIndex.TargetVar} within
the function.}
}
\value{
dataframe of habitat ID and non equalised (R) and equalised (Rg) Phi values
    the R and Rg numbers are both floating point number between -1. and 1.
}
\description{
This function calculates the non and group-equalised versions of the Phi coefficient of
correlation \cite{Tichy and Chytry (2006), De Cáceres and Legendre (2009)}.
}
\details{
It doesn't calculate any p-values by permutation and therefore calculates very quickly
(unlike the permutated version). It is therefore a useful first step and/or useful
when the p-values aren't important.
}
\examples{
    #Builing on the examples from \link[PhiCor]{CorIndex}

    #Using the included dataframe Species1
    Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
    Inputs_species1$nkoverNk
    Species1_AllPhi = CorIndex.all(Inputs_species1)



    #Using the included dataframe Species2
    #If you wanted the analysis to be non-weighted.
    Species2Unweight = Species2[which(Species2$Proportion==1),]
    #This analysis is then the unweighted version of the analysis. SquareID is not needed.
    Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
    names(Inputs_species2)
    Species2_AllPhi = CorIndex.all(Inputs_species2)


}
\references{
Tichy, L., & Chytry, M. (2006). Statistical determination of diagnostic species for site groups of unequal size. Journal of Vegetation Science, 17(6), 809–818.

De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.

Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
}
\seealso{
\link[PhiCor]{PhiCor}
}
\author{
Jordan Chetcuti
}
