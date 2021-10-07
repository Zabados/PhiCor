#' Package for calculation the Weighted Phi Coefficient of association
#'
#' This package can be used to calculate the weighted phi coeficient of
#' association using a weighting for habitats within a location, or
#' can calculate the unweighted original version from
#' \cite{De Cáceres and Legendre (2009)} by using a weighting of one.
#' The package also allows the group-equalised version to be run.
#'
#' The Phi Coefficent of association gives the habitat preference, or
#' site fidelity of habitat sites.
#'
#' @details
#'     The main function in PhiCor is \link[PhiCor]{CorIndex.all.plusP}.
#'     This will calculate all phi corelation of association values and
#'     purmutate to give p-values.
#'
#'     The other Phi coefficent functions are useful in testing or calculating only part
#'     of the full outputs from \link[PhiCor]{CorIndex.all.plusP}:
#'     \itemize{
#'         \item \link[PhiCor]{CorIndex} calculates the inputs that don't change between habitats.
#'         \item \link[PhiCor]{CorIndex.TargetVar} calculates the additional input specific to a habitat.
#'         \item \link[PhiCor]{CorIndex.groupEqual} calculates the group-equalised phi for a habitat.
#'         \item \link[PhiCor]{CorIndex.nonEqual} calculates the non group-equalised phi for a habitat.
#'         \item \link[PhiCor]{CorIndex.xPerm} returns the null phi distribution for a habitat.
#'         \item \link[PhiCor]{CorIndex.all} returns the non and group-equalised phi values for all habitats
#'             without p-values and is therefore quicker.
#'     }
#'     Weighted habitat indicator values (IndVal) can also be calculated, but these have not been tested in
#'     the way that the phi coefficent has been. 
#'     The main function for calculating IndVal is \link[PhiCor]{Indicator.all.plusP}
#'     These functions are useful in testing or calculating only part
#'     of the full outputs from \link[PhiCor]{Indicator.all.plusP}:
#'     \itemize{
#'         \item \link[PhiCor]{Indicator.nonEqual} Calculated the non group-equalised habitat indicator value (IndVal).
#'         \item \link[PhiCor]{Indicator.groupEqual} Calculated the group-equalised species indicator value (IndVal).
#'         \item \link[PhiCor]{Indicator.xPerm} Permutates and calculates for each permutation both the non and 
#'              group-equalised species indicator value (IndVal).
#'         \item 
#'     } 
#'
#'     The functions \link[PhiCor]{lookup} and \link[PhiCor]{FilePathInTree}
#'     might be used elsewhere.
#' @author Jordan Chetcuti
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @keywords internal
#'
"_PACKAGE"
#> [1] "_PACKAGE"






#' Species 1 data of presence absence and weighting for LCM2015 land cover
#'
#' The data is anonymised data created using LCM2015 \cite{Rowland et al., (2017)} and \cite{NBN Atlas website (2018)}
#'
#' @format A data frame with  1,654 rows and 5 variables:
#' @usage The data is useful in examples of how to use the functions.
#' \describe{
#'   \item{LocationID}{A anonymised location}
#'   \item{bhab}{LCM2015 land cover description}
#'   \item{Proportion}{The proportion of a location that is each habitat}
#'   \item{HabId}{An integer ID for each habitat in LCM2015}
#'   \item{Species1}{Binary presence absence for "species1"}
#' }
#' @references Rowland, C. S., Morton, R. D., Carrasco, L., McShane, G., O’Neil, A. W., & Wood, C. M. (2017). Land Cover Map 2015 (vector, GB). NERC Environmental Information Data Centre. NERC Environmental Information Data Centre. doi:/10.5285/6c6c9203-7333-4d96-88ab-78925e7a4e73
#' @references NBN Atlas website. (2018). at https://nbnatlas.org. Accessed 13 April 2018.
"Species1"


#' Species 2 data of presence absence and weighting for LCM2015 land cover
#'
#' The data is anonymised data created using LCM2015 \cite{Rowland et al., (2017)} and \cite{NBN Atlas website (2018)}
#'
#' @format A data frame with  1,005 rows and 5 variables:
#' @usage The data is useful in examples of how to use the functions.
#' \describe{
#'   \item{LocationID}{A anonymised location}
#'   \item{bhab}{LCM2015 land cover description}
#'   \item{Proportion}{The proportion of a location that is each habitat}
#'   \item{HabId}{An integer ID for each habitat in LCM2015}
#'   \item{Species2}{Binary presence absence for "species2"}
#' }
#' @references Rowland, C. S., Morton, R. D., Carrasco, L., McShane, G., O’Neil, A. W., & Wood, C. M. (2017). Land Cover Map 2015 (vector, GB). NERC Environmental Information Data Centre. NERC Environmental Information Data Centre. doi:/10.5285/6c6c9203-7333-4d96-88ab-78925e7a4e73
#' @references NBN Atlas website. (2018). at https://nbnatlas.org. Accessed 13 April 2018.
"Species2"








OutputClass = c(result = data.frame, InVariables = class)
class(OutputClass) = "Output"

print.Output = function(Out){
  print(Out$result)
}

default.Output = function(Out){
  return(Out$result)
}





#' Lookup function for looking up values from one dataframe for another
#'
#' Works very similarly to vlookup in excel.
#' The output is a field in the dataset you are wanting to copy information into.
#' The syntax is something like:
#' MyDataframe$CopiedField = lookup(MyDataframe$joinField, OtherDataframe$joinField, OutherDataframe$DataToCopy)
#'
#' @param MasterDataField vector - a field with common values to those in another dataset you want to
#'     lookup data in.
#' @param AdditionField vector - a field with common values to MasterDataField
#' @param AdditonInputField vector - the data you are importing to the Master dataframe where the
#'     AdditionField matches the MasterDataField
#' @return vector - will reurn the vector corrosponding to the MasterDataField that
#'    are found in the lookup field AdditonInputField
#' @examples
#'    #Lookup DF
#'    HabID = c(1, 2, 3, 4, 5, 6)
#'    Habitat = c("BLeave", "AGrass", "CGrass", "Water", "Bog", "Heath")
#'    Desc = c("Broadleaved and mixed woodland", "Acid grassland","Calcareous grassland", "Inland freshwater", "Sphagnum bog", "Caluna heathland")
#'    OpenClosed = c("Closed", "Open", "Open", "Open","Open","Open")
#'    DomSpecies = c("Oak", "Bent grass", "Festuca ovina", "NA", "Sphagnum Spp.", "Caluna Spp.")
#'    lookupdf = data.frame(HabID, Habitat, Desc, OpenClosed, DomSpecies)
#'
#'    #Feild sites dataframe.
#'    SiteID = c(1,2,3,4,5,6,7,8,9,10)
#'    SiteHabitat = c("BLeave", "AGrass", "BLeave", "AGrass", "Heath", "Bog", "Water", "BLeave", "Water", "Heath")
#'    SitesDF = data.frame(SiteID, SiteHabitat)
#'
#'    #Looking up the descriptions of the habitats.
#'    SitesDF$HabDesc = lookup(SitesDF$SiteHabitat, lookupdf$Habitat, lookupdf$Desc)
#'
#'    #You could lookup other information from the lookup table.
#' @export
lookup = function(MasterDataField, AdditionField, AdditionInputField){
  index <- match(MasterDataField, AdditionField)
  #print(index)
  Field = list()
  count = 1
  for(i in index){
    if(!is.na(i)){
      Field[count] = as.character(AdditionInputField[i])
    }else{
      Field[count] = NA
    }
    count = count + 1
  }
  return(as.character(Field))
}


#' CorIndex is basic initiation of variables that are the same for
#' all habitats in an analysis of multiple habitats in \link[PhiCor]{PhiCor}.
#'
#' This function is usually called by other function, but could be run
#' seperatly to give many of the variables which go into the equations
#' to calculate the phi coefficent. This is useful if manually checking
#' the calculation.
#'
#' @param InDataframe dataframe - containing the data on species presence absence,
#'     habitat and weighting. It also optionally could contain location ID's.
#' @param speciesbinary string - the name of the binary field within InDataframe
#'     that contains binary presence absence for a species at the location.
#' @param weighted string - name of habitat weighting field. For example, proportion
#'     of habitat, but any weighting with all habitats at a location adding to 1.
#' @param group string - ID of habitats or groups of interest.
#' @param SquareID string - optional ID of each location. Useful later for
#'     permutating the locations to calculate p-values.
#' @return a custom class "CorIndexVar" containing all values which are the
#'     same for all habitats and don't need to be recalculated for every
#'     habitat. Including Nk, nk, the habitats "GroupDF", N, K, Ngp, nkoverNk
#'     ng & n. See \cite{De Cáceres and Legendre (2009)}
#'     This is reused to give all the results later.
#' @examples
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'
#' @author Jordan Chetcuti
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
CorIndex <- function(InDataframe, speciesbinary,weighted, group, SquareID = NULL){
  speciesbinary = InDataframe[,speciesbinary]
  weighted = InDataframe[,weighted]
  group = InDataframe[,group]
  NextMethod("CorIndex")
}

#' The default method of \link[PhiCor]{CorIndex}.
#'
#' This should not be used independently. See \link[PhiCor]{CorIndex} for usage.
#' @usage NULL
#' @export
CorIndex.default <- function(InDataframe, speciesbinary,weighted, group, SquareID = NULL){
  rtrn <- list()
  class(rtrn) <- "CorIndexVar"
  ### Weighting the variables Np and np ( calculated for every habitat as Nk and nk and then
  ###later individually selected as Np and np)
  Nk = aggregate(weighted ~ group, InDataframe, sum)
  Nk = setNames(Nk, c("group", "Nk"))
  nk = aggregate(speciesbinary * weighted ~ group, InDataframe, sum)
  nk = setNames(nk, c("group", "nk"))
  rtrn$Nk = Nk
  rtrn$nk = nk
  GroupDF = Nk
  GroupDF[,2] =NULL
  rtrn$GroupDF = GroupDF
  N = sum(weighted)
  rtrn$N = N
  K = length(Nk$group)
  rtrn$K = K
  Ngp = N/K
  rtrn$Ngp = Ngp
  nkoverNk = nk/Nk
  nkoverNk = setNames(nkoverNk, c("group", "nkoverNk"))
  rtrn$nkoverNk = nkoverNk
  sumnkoverNk = sum(nkoverNk$nkoverNk)
  rtrn$sumnkoverNk = sumnkoverNk
  ng = Ngp*sumnkoverNk
  rtrn$ng = ng
  n = sum(nk$nk)
  rtrn$n = n
  return(rtrn)
}


#' CorIndex.TargetVar is creates the variables that are particular to a
#' habitat that were not created in the \link[PhiCor]{CorIndex}
#'
#' This function is usually called by other function, but could be run
#' seperatly to give many of the variables which go into the equations
#' to calculate the phi coefficent. This is useful if manually checking
#' the calculation.
#'
#' @param CorIndexVarInput is the custom class "CorIndexVar" output from \link[PhiCor]{CorIndex}
#' @param targetgroup Integer - ID of the habitat of interest.
#' @return a custom class "CorIndexVar" containing all values which are the
#'     same for all habitats and don't need to be recalculated for every
#'     habitat and the ones for the specific target habitat. Including those from
#'     CorIndex Nk, nk, the habitats "GroupDF", N, K, Ngp, nkoverNk ng & n and
#'     aditionally np, Np & ngp  See \cite{De Cáceres and Legendre (2009)}
#'     This is reused to give all the results later.
#'
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of CorIndex. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = CorIndex.TargetVar(Inputs_species1, 1)
#'
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'
#'     Looking at habitat 4
#'     Inputs_species2_hab4 = CorIndex.TargetVar(Inputs_species2, 4)
#'
#'
#' @seealso \link[PhiCor]{CorIndex}
#' @author Jordan Chetcuti
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
CorIndex.TargetVar <- function(CorIndexVarInput, targetgroup){
  rtrn <- list()
  class(rtrn) <- "CorIndexVar"
  rtrn = CorIndexVarInput
  nk = CorIndexVarInput$nk
  Nk = CorIndexVarInput$Nk
  Ngp = CorIndexVarInput$Ngp
  np = nk$nk[nk$group==targetgroup]
  Np = Nk$Nk[Nk$group==targetgroup]
  rtrn$np = np
  rtrn$Np = Np
  ngp = Ngp*np/Np
  rtrn$ngp = ngp
  return(rtrn)
}

#' Calculated the group-equalised Phi Coefficient of correlation.
#'
#' This function calculates the group-equalised version of the Phi coefficient of
#' correlation \cite{Tichy and Chytry (2006), De Cáceres and Legendre (2009)}. This is
#' useful alone to calculate the group-equalised version for a single habitat of interest
#' among multiple other habitats.
#'
#' @param CorIndexTargetVarInput a custom class "CorIndexVar" containing all values from
#'     from both \link[PhiCor]{CorIndex} and \link[PhiCor]{CorIndex.TargetVar}
#' @return Phi - a floating point number between -1. and 1.
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex.TargetVar}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of CorIndex. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = CorIndex.TargetVar(Inputs_species1, 1)
#'
#'     species1_hab1_GE_Phi = CorIndex.groupEqual(Inputs_species1_hab1)
#'
#'
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'
#'     #Looking at habitat 4
#'     Inputs_species2_hab4 = CorIndex.TargetVar(Inputs_species2, 4)
#'
#'     species2_hab4_GE_Phi = CorIndex.groupEqual(Inputs_species2_hab4)
#'
#'
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @references
#'     Tichy, L., & Chytry, M. (2006). Statistical determination of diagnostic species for site groups of unequal size. Journal of Vegetation Science, 17(6), 809–818.
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
CorIndex.groupEqual <- function(CorIndexTargetVarInput){
  N = CorIndexTargetVarInput$N
  ngp = CorIndexTargetVarInput$ngp
  ng = CorIndexTargetVarInput$ng
  Ngp = CorIndexTargetVarInput$Ngp
  Requ = (N*ngp - ng*Ngp)/sqrt((N*ng - ng^2)*(N*Ngp - Ngp^2))
  return(Requ)
}

#' Calculated the non group-equalised Phi Coefficient of correlation.
#'
#' This function calculates the non group-equalised version of the Phi coefficient of
#' correlation \cite{De Cáceres and Legendre (2009)}. This is
#' useful alone to calculate the non group-equalised version for a single habitat of interest
#' among multiple other habitats.
#'
#' @param CorIndexTargetVarInput a custom class "CorIndexVar" containing all values from
#'     from both \link[PhiCor]{CorIndex} and \link[PhiCor]{CorIndex.TargetVar}
#' @return Phi - a floating point number between -1. and 1.
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex.TargetVar}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of CorIndex. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = CorIndex.TargetVar(Inputs_species1, 1)
#'
#'     species1_hab1_Phi = CorIndex.nonEqual(Inputs_species1_hab1)
#'
#'
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'
#'     #Looking at habitat 4
#'     Inputs_species2_hab4 = CorIndex.TargetVar(Inputs_species2, 4)
#'
#'     species2_hab4_Phi = CorIndex.nonEqual(Inputs_species2_hab4)
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#'     Tichy, L., & Chytry, M. (2006). Statistical determination of diagnostic species for site groups of unequal size. Journal of Vegetation Science, 17(6), 809–818.
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
CorIndex.nonEqual <- function(CorIndexTargetVarInput){
  N = CorIndexTargetVarInput$N
  np = CorIndexTargetVarInput$np
  n = CorIndexTargetVarInput$n
  Np = CorIndexTargetVarInput$Np
  Rnon = (N*np - n*Np)/sqrt((N*n - n^2)*(N*Np - Np^2))
  return(Rnon)
}

#' Calculated both the non and group-equalised Phi Coefficient of correlation for all habitats.
#'
#' This function calculates the non and group-equalised versions of the Phi coefficient of
#' correlation \cite{Tichy and Chytry (2006), De Cáceres and Legendre (2009)}.
#'
#' It doesn't calculate any p-values by permutation and therefore calculates very quickly
#' (unlike the permutated version). It is therefore a useful first step and/or useful
#' when the p-values aren't important.
#'
#' @param CorIndexVarInput a custom class "CorIndexVar" containing all values from
#'     from only \link[PhiCor]{CorIndex}. It calls \link[PhiCor]{CorIndex.TargetVar} within
#'     the function.
#' @return dataframe of habitat ID and non equalised (R) and equalised (Rg) Phi values
#'     the R and Rg numbers are both floating point number between -1. and 1.
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex}
#'
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     Species1_AllPhi = CorIndex.all(Inputs_species1)
#'
#'
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'     Species2_AllPhi = CorIndex.all(Inputs_species2)
#'
#'
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @references
#'     Tichy, L., & Chytry, M. (2006). Statistical determination of diagnostic species for site groups of unequal size. Journal of Vegetation Science, 17(6), 809–818.
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
CorIndex.all <- function(CorIndexVarInput){
  OutDF = CorIndexVarInput$GroupDF
  OutDF$R = 0.0000
  OutDF$Rg = 0.0000
  for(i in 1:nrow(OutDF)) {
    targetVariables = CorIndex.TargetVar(CorIndexVarInput,OutDF[i,1])
    print(OutDF[i,1])
    OutDF[i,2]= CorIndex.nonEqual(targetVariables)
    OutDF[i,3]= CorIndex.groupEqual(targetVariables)
  }
  return(OutDF)
}

CorIndex.APerm <- function(InDataframe, FieldToPerm, SquareID = NULL){
  copyDF = InDataframe
  if(is.null(SquareID)){
    copyDF[,FieldToPerm] = sample(copyDF[,FieldToPerm], nrow(copyDF), replace = FALSE, prob = NULL)
  }else{
    AggCopyDF = aggregate(InDataframe[,FieldToPerm] ~ InDataframe[,SquareID], InDataframe, min)
    AggCopyDF = setNames(AggCopyDF, c(SquareID, FieldToPerm))
    AggCopyDF[,FieldToPerm] = sample(AggCopyDF[,FieldToPerm], nrow(AggCopyDF), replace = FALSE, prob = NULL)
    copyDF[,FieldToPerm] = as.numeric(lookup(copyDF[,SquareID], AggCopyDF[,SquareID], AggCopyDF[,2]))
  }
  return(copyDF)
}


CorIndex.toroidalAPerm <- function(InDataframe, FieldToPerm, SquareID = NULL){
  copyDF = InDataframe
  
  
  if(is.null(SquareID)){
    rows = nrow(copyDF)
    MoveOn = sample(0:rows, 1)
    
    copyDF$Index = row.names(copyDF)
    copyDF$Index = as.numeric(copyDF$Index) + MoveOn
    
    copyDF$Index = ifelse(copyDF$Index>rows,copyDF$Index - rows, copyDF$Index )
    
    temp = copyDF[,c(FieldToPerm, "Index")]
    
    temp = temp[order(temp$Index),]
    
    copyDF[,FieldToPerm] =  temp[,FieldToPerm]
    copyDF$Index = NULL
    
    #copyDF[,FieldToPerm] = sample(copyDF[,FieldToPerm], nrow(copyDF), replace = FALSE, prob = NULL)
  }else{

    AggCopyDF = aggregate(InDataframe[,FieldToPerm] ~ InDataframe[,SquareID], InDataframe, min)
    AggCopyDF = setNames(AggCopyDF, c(SquareID, FieldToPerm))
    
    
    
    
    rows = nrow(AggCopyDF)
    MoveOn = sample(0:rows, 1)
    
    AggCopyDF$Index = row.names(AggCopyDF)
    AggCopyDF$Index = as.numeric(AggCopyDF$Index) + MoveOn
    
    AggCopyDF$Index = ifelse(AggCopyDF$Index>rows,AggCopyDF$Index - rows, AggCopyDF$Index )
    
    temp = AggCopyDF[,c(FieldToPerm, "Index")]
    
    temp = temp[order(temp$Index),]
    
    AggCopyDF[,FieldToPerm] =  temp[,FieldToPerm]
    AggCopyDF$Index = NULL
    
    
    
    
    
    
    
    
    
    
    AggCopyDF[,FieldToPerm] = sample(AggCopyDF[,FieldToPerm], nrow(AggCopyDF), replace = FALSE, prob = NULL)
    copyDF[,FieldToPerm] = as.numeric(lookup(copyDF[,SquareID], AggCopyDF[,SquareID], AggCopyDF[,2]))
  }
  return(copyDF)
}


#' Permutates and calculates for each permutation both the non and group-equalised
#' Phi coefficent.
#'
#' This function purmutaes and calculates the non and group-equalised versions of the Phi coefficient of
#' correlation \cite{Tichy and Chytry (2006), De Cáceres and Legendre (2009)}.
#' It is useful as it gives the distribution of phi values if the species are distributed
#' randomly across the locations, the phi value for the actual data is compared to this
#' distribution to calculate the p-value. The distribution is a useful output to visualise
#' the distribution.
#'
#' @param InDataframe dataframe - containing the data on species presence absence,
#'     habitat and weighting. It also optionally could contain location ID's.
#' @param FieldToPerm string - the name of the binary field within InDataframe
#'     that contains binary presence absence for a species at the location.
#' @param weighted string - name of habitat weighting field. For example, proportion
#'     of habitat, but any weighting with all habitats at a location adding to 1.
#' @param groupsField string - ID of habitats or groups of interest.
#' @param targetgroup Integer - ID of the habitat of interest.
#' @param x Integer - optional - the number of permutation, set as 1000 by default.
#' @param SquareID string - optional if non-weighted, but MUST be used if a weighting
#'     is used. ID of each location. Useful later for permutating the locations to
#'     calculate p-values.
#' @param toroidal Boolean - If you suspect that your data is spatially auto-correlated you should use 
#'     a toroidal permutation.
#' @return Two data frames, one for the group-equalised $`GroupEqualPermList` and one
#'     for the non group-equalised $NonGroupEqualPermList. It also returns the number of
#'     permutations $Permutations.
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex.nonEqual}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of CorIndex. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = CorIndex.TargetVar(Inputs_species1, 1)
#'
#'     species1_hab1_Phi = CorIndex.nonEqual(Inputs_species1_hab1)
#'
#'     #histogram showing the output of CorIndex.xPerm for habitat 1 in comparison to the phi value from actual data.
#'     hist(Hab1_PermutatedPhi$NonGroupEqualPermList, main = "Non group-equalised distribution", xlab = "Phi values", breaks = 100)
#'     abline(v=species1_hab1_Phi, col="red")
#'     text((species1_hab1_Phi + 0.006), 15, "Calulated phi value", col = "red", srt=90)
#'
#'
#'
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @export
CorIndex.xPerm <- function(InDataframe, FieldToPerm, weighted, groupsField, targetGroup, x = 1000, SquareID = NULL, toroidal = FALSE){
  GroupEqualPermList = c()
  NonGroupEqualPermList = c()
  for(i in 1:x){
    
    if (toroidal){
      Perm1 = CorIndex.toroidalAPerm(InDataframe, FieldToPerm, SquareID)
    }else{
      Perm1 = CorIndex.APerm(InDataframe, FieldToPerm, SquareID)
    }
    
    
    UnchangePerm = CorIndex(Perm1, FieldToPerm, weighted, groupsField)
    VarPerm1 = CorIndex.TargetVar(UnchangePerm, targetGroup)
    GroupEqualPermList = c(GroupEqualPermList, CorIndex.groupEqual(VarPerm1))
    NonGroupEqualPermList = c(NonGroupEqualPermList, CorIndex.nonEqual(VarPerm1))
  }
  out = list("GroupEqualPermList" = GroupEqualPermList, "NonGroupEqualPermList" = NonGroupEqualPermList, "Permutations" = x)
  return(out)
}

#' Calculated both the non and group-equalised Phi Coefficient of correlation for all habitats
#' and the p-values.
#'
#' This function calculates the non and group-equalised versions of the Phi coefficient of
#' correlation \cite{Tichy and Chytry (2006), De Cáceres and Legendre (2009)}.
#'
#' It calculate p-values by permutation and therefore calculates take a bit of time.
#' This is the main way of useing the PhiCor package and give full outputs.
#'
#' @param InDataframe dataframe - containing the data on species presence absence,
#'     habitat and weighting. It also optionally could contain location ID's.
#' @param speciesbinary string - the name of the binary field within InDataframe
#'     that contains binary presence absence for a species at the location.
#' @param weighted string - name of habitat weighting field. For example, proportion
#'     of habitat, but any weighting with all habitats at a location adding to 1.
#' @param group string - ID of habitats or groups of interest.
#' @param numberIteration Integer - optional - the number of permutation, set as 1000 by default.
#' @param SquareID string - optional if non-weighted, but MUST be used if a weighting
#'     is used. ID of each location. Useful later for permutating the locations to
#'     calculate p-values.
#' @param toroidal Boolean - If you suspect that your data is spatially auto-correlated you should use 
#'     a toroidal permutation.
#' @return a custom class "CorIndexVar" containing all values which are the
#'     same for all habitats and don't need to be recalculated for every
#'     habitat. Including Nk, nk, the habitats "GroupDF", N, K, Ngp, nkoverNk
#'     ng & n. See \cite{De Cáceres and Legendre (2009)}. To this is added another custom
#'     class "Output" as $`result` which gives the habitat ID, R (non group-equalised phi)
#'     pR (the p-value for R), Rg (group-equalised phi) and pRg (the p-value for Rg). This
#'     "Output" class is the main result, the CorIndexVar is included as metadata of the
#'     input values.
#' @examples
#'     #Using the included dataframe Species1
#'     Species1_AllPhiPvalues = CorIndex.all.plusP(InDataframe = Species1, speciesbinary = "Species1", weighted = 'Proportion', group = 'HabId', SquareID = 'LocationID')
#'     names(Species1_AllPhiPvalues)
#'     Species1_AllPhiPvalues$`result`
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @references
#'     Tichy, L., & Chytry, M. (2006). Statistical determination of diagnostic species for site groups of unequal size. Journal of Vegetation Science, 17(6), 809–818.
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
CorIndex.all.plusP <- function(InDataframe, speciesbinary, weighted, group, numberIteration = 1000, SquareID = NULL, toroidal = FALSE){
  CorIndexVarInput = CorIndex(InDataframe, speciesbinary, weighted, group, SquareID)
  OutDF = CorIndexVarInput$GroupDF
  OutDF$R = 0.0000000000
  OutDF$pR = 0.0000000000
  OutDF$Rg = 0.0000000000
  OutDF$pRG = 0.0000000000
  for(i in 1:nrow(OutDF)) {
    targetVariables = CorIndex.TargetVar(CorIndexVarInput,OutDF[i,1])
    print(OutDF[i,1])
    Distributions = CorIndex.xPerm(InDataframe, speciesbinary,weighted,group,OutDF[i,"group"], numberIteration, SquareID, toroidal)
    OutDF[i,2]= CorIndex.nonEqual(targetVariables)
    #Calc each of the probabilities.
    length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList >= OutDF[i,2]])/numberIteration
    length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList <= OutDF[i,2]])/numberIteration
    OutDF[i,3] = min(length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList >= OutDF[i,2]])/numberIteration , length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList <= OutDF[i,2]])/numberIteration)
    OutDF[i,4]= CorIndex.groupEqual(targetVariables)
    OutDF[i,5] = min(length(Distributions$GroupEqualPermList[Distributions$GroupEqualPermList >= OutDF[i,4]])/numberIteration , length(Distributions$GroupEqualPermList[Distributions$GroupEqualPermList <= OutDF[i,4]])/numberIteration)
    OutputClass$result = OutDF
    OutputClass$InVariables = CorIndexVarInput
  }
  return(OutputClass)
}



#' Calculated the non group-equalised habitat indicator value (IndVal).
#'
#' This function calculates the non group-equalised version of the habitat indicator value (IndVal)
#' \cite{De Cáceres and Legendre (2009)}. This is
#' useful alone to calculate the non group-equalised version for a single habitat of interest
#' among multiple other habitats.
#'
#' @param CorIndexTargetVarInput a custom class "CorIndexVar" containing all values from
#'     from both \link[PhiCor]{CorIndex} and \link[PhiCor]{CorIndex.TargetVar}
#' @return IndVal - a floating point number between 0 and 1.
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex.TargetVar}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of Indicator. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = Indicator.TargetVar(Inputs_species1, 1)
#'
#'     species1_hab1_Phi = Indicator.nonEqual(Inputs_species1_hab1)
#'
#'
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'
#'     #Looking at habitat 4
#'     Inputs_species2_hab4 = Indicator.TargetVar(Inputs_species2, 4)
#'
#'     species2_hab4_IndVal = Indicator.nonEqual(Inputs_species2_hab4)
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
Indicator.nonEqual <- function(CorIndexTargetVarInput){
  N = CorIndexTargetVarInput$N
  np = CorIndexTargetVarInput$np
  n = CorIndexTargetVarInput$n
  Np = CorIndexTargetVarInput$Np
  ########
  Rnon = sqrt((np/n)*(np/Np))
  return(Rnon)
}

#' Calculated the group-equalised species indicator value (IndVal).
#'
#' This function calculates the group-equalised version of the species indicator value (IndVal)
#' \cite{De Cáceres and Legendre (2009)}. This is
#' useful alone to calculate the group-equalised version for a single habitat of interest
#' among multiple other habitats.
#'
#' @param CorIndexTargetVarInput a custom class "CorIndexVar" containing all values from
#'     from both \link[PhiCor]{CorIndex} and \link[PhiCor]{CorIndex.TargetVar}
#' @return IndVal - a floating point number between 0 and 1.
#' @examples
#'     #Builing on the examples from \link[PhiCor]{CorIndex.TargetVar}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of Indicator. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = CorIndex.TargetVar(Inputs_species1, 1)
#'
#'     species1_hab1_GE_IndVal = Indicator.groupEqual(Inputs_species1_hab1)
#'
#'
#'
#'     #Using the included dataframe Species2
#'     #If you wanted the analysis to be non-weighted.
#'     Species2Unweight = Species2[which(Species2$Proportion==1),]
#'     #This analysis is then the unweighted version of the analysis. SquareID is not needed.
#'     Inputs_species2 =  CorIndex(InDataframe = Species2Unweight, speciesbinary = 'Species2', weighted = 'Proportion', group = 'HabId')
#'     names(Inputs_species2)
#'
#'     #Looking at habitat 4
#'     Inputs_species2_hab4 = Indicator.TargetVar(Inputs_species2, 4)
#'
#'     species2_hab4_GE_IndVal = Indicator.groupEqual(Inputs_species2_hab4)
#'
#'
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
Indicator.groupEqual <- function(CorIndexTargetVarInput){
  N = CorIndexTargetVarInput$N
  ngp = CorIndexTargetVarInput$ngp
  ng = CorIndexTargetVarInput$ng
  Ngp = CorIndexTargetVarInput$Ngp
  
  
  SumNks = sum(CorIndexTargetVarInput$nk/CorIndexTargetVarInput$Nk)
  
  ################
  Requ = sqrt(((ngp/Ngp)/SumNks)*(ngp/Ngp))
  return(Requ)
}





#' Permutates and calculates for each permutation both the non and group-equalised
#' species indicator value (IndVal).
#'
#' This function purmutaes and calculates the non and group-equalised versions of species indicator 
#' value (IndVal) \cite{De Cáceres and Legendre (2009)}.
#' It is useful as it gives the distribution of IndVal values if the species are distributed
#' randomly across the locations, the IndVal value for the actual data is compared to this
#' distribution to calculate the p-value. The distribution is a useful output to visualise
#' the distribution.
#'
#' @param InDataframe dataframe - containing the data on species presence absence,
#'     habitat and weighting. It also optionally could contain location ID's.
#' @param FieldToPerm string - the name of the binary field within InDataframe
#'     that contains binary presence absence for a species at the location.
#' @param weighted string - name of habitat weighting field. For example, proportion
#'     of habitat, but any weighting with all habitats at a location adding to 1.
#' @param groupsField string - ID of habitats or groups of interest.
#' @param targetgroup Integer - ID of the habitat of interest.
#' @param x Integer - optional - the number of permutation, set as 1000 by default.
#' @param SquareID string - optional if non-weighted, but MUST be used if a weighting
#'     is used. ID of each location. Useful later for permutating the locations to
#'     calculate p-values.
#' @param toroidal Boolean - If you suspect that your data is spatially auto-correlated you should use 
#'     a toroidal permutation.
#' @return Two data frames, one for the group-equalised $`GroupEqualPermList` and one
#'     for the non group-equalised $NonGroupEqualPermList. It also returns the number of
#'     permutations $Permutations.
#' @examples
#'     #Building on the examples from \link[PhiCor]{Indicator.nonEqual}
#'     #Using the included dataframe Species1
#'     Inputs_species1 = CorIndex(InDataframe = Species1, speciesbinary = 'Species1', weighted = 'Proportion', group = 'HabId', 'LocationID')
#'     Inputs_species1$nkoverNk
#'     #Using CorIndex.TargetVar, the input is the output of Indicator. In this case we are interested in
#'     #habitat 1 (using the ID for the habitat)
#'     Inputs_species1_hab1 = CorIndex.TargetVar(Inputs_species1, 1)
#'
#'     species1_hab1_Indval = Indicator.nonEqual(Inputs_species1_hab1)
#'
#'     #histogram showing the output of Indicator.xPerm for habitat 1 in comparison to the IndVal value from actual data.
#'     hist(Hab1_PermutatedPhi$NonGroupEqualPermList, main = "Non group-equalised distribution", xlab = "IndVal values", breaks = 100)
#'     abline(v=species1_hab1_IndVal, col="red")
#'     text((species1_hab1_IndVal + 0.006), 15, "Calulated IndVal value", col = "red", srt=90)
#'
#'
#'
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @export
Indicator.xPerm <- function(InDataframe, FieldToPerm, weighted, groupsField, targetGroup, x = 1000, SquareID = NULL, toroidal = FALSE){
  GroupEqualPermList = c()
  NonGroupEqualPermList = c()
  for(i in 1:x){
    if (toroidal){
      Perm1 = CorIndex.toroidalAPerm(InDataframe, FieldToPerm, SquareID)
    }else{
      Perm1 = CorIndex.APerm(InDataframe, FieldToPerm, SquareID)
    }
    
    
    UnchangePerm = CorIndex(Perm1, FieldToPerm, weighted, groupsField)
    VarPerm1 = CorIndex.TargetVar(UnchangePerm, targetGroup)
    
    #########
    GroupEqualPermList = c(GroupEqualPermList, Indicator.groupEqual(VarPerm1))
    NonGroupEqualPermList = c(NonGroupEqualPermList, Indicator.nonEqual(VarPerm1))
  }
  out = list("GroupEqualPermList" = GroupEqualPermList, "NonGroupEqualPermList" = NonGroupEqualPermList, "Permutations" = x)
  return(out)
}


#' Calculated both the non and group-equalised IndVal for all habitats
#' and the p-values.
#'
#' This function calculates the non and group-equalised versions of an habitat indicator value (IndVal)
#' correlation \cite{De Cáceres and Legendre (2009)}.
#'
#' It calculate p-values by permutation and therefore calculates take a bit of time.
#'
#' @param InDataframe dataframe - containing the data on species presence absence,
#'     habitat and weighting. It also optionally could contain location ID's.
#' @param speciesbinary string - the name of the binary field within InDataframe
#'     that contains binary presence absence for a species at the location.
#' @param weighted string - name of habitat weighting field. For example, proportion
#'     of habitat, but any weighting with all habitats at a location adding to 1.
#' @param group string - ID of habitats or groups of interest.
#' @param numberIteration Integer - optional - the number of permutation, set as 1000 by default.
#' @param SquareID string - optional if non-weighted, but MUST be used if a weighting
#'     is used. ID of each location. Useful later for permutating the locations to
#'     calculate p-values.
#' @param toroidal Boolean - If you suspect that your data is spatially auto-correlated you should use 
#'     a toroidal permutation.
#' @return a custom class "CorIndexVar" containing all values which are the
#'     same for all habitats and don't need to be recalculated for every
#'     habitat. Including Nk, nk, the habitats "GroupDF", N, K, Ngp, nkoverNk
#'     ng & n. See \cite{De Cáceres and Legendre (2009)}. To this is added another custom
#'     class "Output" as $`result` which gives the habitat ID, I (non group-equalised IndVal)
#'     pI (the p-value for I), Ig (group-equalised IndVal) and Ig (the p-value for Ig). This
#'     "Output" class is the main result, the CorIndexVar is included as metadata of the
#'     input values.
#' @examples
#'     #Using the included dataframe Species1
#'     Species1_AllIndValvalues = Indicator.all.plusP(InDataframe = Species1, speciesbinary = "Species1", weighted = 'Proportion', group = 'HabId', SquareID = 'LocationID')
#'     names(Species1_AllIndValvalues)
#'     Species1_AllIndValvalues$`result`
#' @seealso \link[PhiCor]{PhiCor}
#' @author Jordan Chetcuti
#' @references
#'     De Cáceres, M., & Legendre, P. (2009). Associations between species and groups of sites: indices and statistical inference. Ecology, 90(12), 3566–3574.
#' @references
#'     Chetcuti, J., Kunin, W. E. & Bullock, J. M. A weighting method to improve habitat association analysis: tested on British carabids. Ecography (Cop.). 42, 1395–1404 (2019).
#' @export
Indicator.all.plusP <- function(InDataframe, speciesbinary, weighted, group, numberIteration = 1000, SquareID = NULL, toroidal = NULL){
  CorIndexVarInput = CorIndex(InDataframe, speciesbinary, weighted, group, SquareID)
  OutDF = CorIndexVarInput$GroupDF
  OutDF$I = 0.0000000000
  OutDF$pI = 0.0000000000
  OutDF$Ig = 0.0000000000
  OutDF$pIG = 0.0000000000
  
  
  for(i in 1:nrow(OutDF)) {
    targetVariables = CorIndex.TargetVar(CorIndexVarInput,OutDF[i,1])
    
    
    print(OutDF[i,1])
    
    Distributions = Indicator.xPerm(InDataframe, speciesbinary,weighted,group,i, numberIteration, SquareID, toroidal)
    
    
    OutDF[i,2]= Indicator.nonEqual(targetVariables)
    
    
    
    #Calc each of the probabilities.
    length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList >= OutDF[i,2]])/numberIteration
    length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList <= OutDF[i,2]])/numberIteration
    OutDF[i,3] = min(length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList >= OutDF[i,2]])/numberIteration , length(Distributions$NonGroupEqualPermList[Distributions$NonGroupEqualPermList <= OutDF[i,2]])/numberIteration)
    OutDF[i,4]= Indicator.groupEqual(targetVariables)
    OutDF[i,5] = min(length(Distributions$GroupEqualPermList[Distributions$GroupEqualPermList >= OutDF[i,4]])/numberIteration , length(Distributions$GroupEqualPermList[Distributions$GroupEqualPermList <= OutDF[i,4]])/numberIteration)
    
    OutputClass$result = OutDF
    OutputClass$InVariables = CorIndexVarInput
  }
  return(OutputClass)
}


#' Looks for a filename in the current folder and then expands up and along the file tree
#' looking for the filename.
#'
#' This isn't really part of calculating the Phi coefficient, but is useful in finding a
#' file which has been reorganised within a project. Imperfect as will return multiple
#' paths if the filename is used twice locally.
#'
#' @param StringFileName string filename
#' @return Filepath including the input StringFileName for all files found in the vicinity.
#'     If it finds multiple files it returns a list of paths which can be accessed individually
#'     using [1], [2]...etc.
#' @examples
#'     #Will return as "character" all R related files in the directory you are currently working on
#'     PathRfiles = FilePathInTree(".R")
#'     #If there is only one file PathRfiles would return it, but with multiple PathRfiles[1] would return the first or only.
#'     PathRfiles[1]
#'
#' @export
FilePathInTree = function(StringFileName){
  OldPath = getwd()
  while(length(dir(getwd(), StringFileName,recursive=TRUE))==0){
    setwd('..')
  }
  FilePath = file.path(getwd(),dir(getwd(), StringFileName,recursive=TRUE))
  setwd(OldPath)
  rm(OldPath)
  return(FilePath)
}
