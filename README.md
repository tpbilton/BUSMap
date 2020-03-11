# BUSMap

[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

Bayesian genotyping Uncertainty with Sequencing data and linkage MAPping (BUSMap).

An R package for constructing linkage maps in full-sib populations using a Bayesian hierarchical hidden Markov model (HMM).

### Installation:

The easiest way to install BUSMap in R is using the devtools package.

```
install.packages("devtools")
devtools::install_github("tpbilton/BUSMap")
```

### Example:

We give a simple example to illustrate how to use BUSMap. We fit the Bayesian model to the manuka data found in Bilton et al. (2018) using the combined low depth and high depth SNP sets.

The data can be loaded using the following code:

```
data(manuka)         # load data in BUSMap   
ref <- manuka$ref    # matrix of reference read counts 
alt <- manuka$alt    # matrix of alternate read counts
parHap <- parHap     # Matrix of parental haplotypes
```

A linkage map can be computed via the Bayesian hierarchical HMM using the `computeMapSeq` function:
```
## Simulation parameters
burn = 5000  # burn-in period
iter = 25000 # number of iterations (excluding burn-in)  
nchain = 3   # number of chains

## Construct the linkage map using Bayeisan hierarchical HMM
est <- computeMapSeq(ref, alt, parHap, iter=iter, burnin=burn, chains=nchain)
str(est)
```
The output is a list of matrices of the posterior samples. Each list is samples from one chain, the rows of each matrix represent the iterations in the MH algorithm and the columns represent the parameters.  

### Citation:
Bilton, T.P., Schofield, M.R., Black, M.A., Chagn&#233;, D., Wilcox, P.L., & Dodds, K.G. (2018). Accounting for errors in low coverage high-throughput sequencing data when constructing genetic maps using biparental outcrossed populations. *Genetics*, *209*(1), 65--76. doi:[10.1534/genetics.117.300627](http://www.genetics.org/content/209/1/65) 

### Funding:
The initial development of this package was partially funded by the Ministry of Business, Innovation and Employment via its funding of the “Genomics for Production & Security in a Biological Economy” programme (Contract ID C10X1306).


