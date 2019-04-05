# PhiCor

----------------
A package for calculating the phi Coefficient of Association as a means of determining habitat association/preference.
You can also use the package to calculate indicator values (IndVal).

### Cite
----------------
Chetcuti, J. (2018). PhiCor: PhiCor: Calculating the Phi coefficient of Association. R package version 0.1.0. https://doi.org/10.5281/zenodo.1407053

### References
----------------
Chetcuti, J. et al. 2019. A weighting method to improve habitat association analysis: tested on British carabids. - Ecography (Cop.).: ecog.04295. https://doi.org/10.1111/ecog.04295

### Data
----------------
Chetcuti, J.; Kunin, W.E.; Bullock, J.M. (2019). Ground beetle (Carabidae) associations with UK Land Cover Map habitats. NERC Environmental Information Data Centre. https://doi.org/10.5285/ce0a6690-9277-4880-a20a-b30477bf8646

Example results: https://shiny-apps.ceh.ac.uk/CarabidData/

### Installation
----------------

Use the `devtools` package to **install**.

```r
# install.packages("devtools")
# NOTE: If you have not installed devtools before you will need to restart your R session

library(devtools)

# Install PhiCor
install_github('Zabados/PhiCor')

# Load PhiCor
library(PhiCor)
```
