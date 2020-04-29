##########################################################################
# Bayesian genotyping Uncertainty with Sequencing data and linkage MAPping (BUSMap)
# Copyright 2019-2020 Timothy P. Bilton <timothy.bilton@agresearch.co.nz>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#########################################################################
#' Construct a linkage map using high-throughput sequencing data
#' 
#' Function that implements a Metropolis-Hastings (MH) algorithm for obtaining posterior samples 
#' of the Bayesian hierarchical hidden Markov model (HMM) for linkage maps in full-sib families 
#' with high-throughput sequencing data. Note that this function only constructs a linkage map for a single chromsome
#' in a single full-sib family.
#' 
#' The matrices for \code{ref} and \code{alt} must have the samples (or inidividuals) indexed by the rows
#' and the SNPs by the columns. Since this function is for full-sib families, only the read counts for the 
#' offspring should be included (the parental information is contained in \code{parhap}).
#' 
#' The \code{parhap} argument must be a matrix (4 rows and \eqn{M} columns, where \eqn{M} is the number of SNPs) containing the haplotypes for the parents. The entries 
#' of this matrix can be anything but there must be excatly two unique entires, one to denote the reference allele and
#' one to denote the alternate allele. The first two rows specify the haplotypes for the maternal parent and the third and
#' fourth rows specify the haplotypes for the paternal parent. See examples for one way to specify this matrix.
#' 
#' The \code{seed} argument allows for reproducibility and specifies the starting seed to when simulating random
#' values in the MCMC chains. 
#' 
#' The \code{initval} allows initial values for each chain to be specified. This must a list of vectors (each with length \eqn{2M+3}) with the number of elements
#' in the list being equal to the number of chains specified by \code{chains}. This argument should only be used by experienced users.  
#' 
#' @usage
#' computeMapSeq(ref, alt, parHap, iter=30000, burnin=5000, chains=3, seed=1, cores=chains, initval=NULL)
#' 
#' @param ref Non-negative integer matrix containing the read counts for the reference allele.
#' @param alt Non-negative integer matrix containing the read counts for the alternate allele.
#' @param parhap Matrix of parental haplotypes.
#' @param iter Integer value of the number of iterations in the MH algorithm (excluding burn-in period) for each chain.
#' @param burnin Integer value of the number of iterations in the burn-in/adaptive phase of the MH algorithm for each chain.
#' @param chains Integer value giving the number of parallel chains to use.  
#' @param seed Numeric value giving the seed.
#' @param cores Integer value giving the number of cores for parallelization of the MH chains.
#' @param initval List giving the starting values for the MH algorithm value.
#' 
#' @name computeMapSeq
#' @author Timothy P. Bilton
#' @examples
#' 
#' #### Load Manuka data from Bilton et al., (2018)
#' data(manuka)
#' ref <- manuka$ref       # matrix of reference allele counts
#' alt <- manuka$alt       # matrix of alternate allele counts
#' parHap <- manuka$parHap # Matrix of parental haplotypes
#' 
#' computeMapSeq(ref, alt, parHap)
#'
#' @export computeMapSeq
computeMapSeq <- function(ref, alt, parhap, iter=30000, burnin=5000, chains=3, seed=as.numeric(Sys.time()), cores=chains, initval=NULL){
  ## Do some checks
  if(!is.matrix(ref) || !is.numeric(ref) || any(is.na(ref)) || any(ref != round(ref)) || any(ref < 0))
    stop("Argument `ref` is invalid")
  if(!is.matrix(alt) || !is.numeric(alt) || any(is.na(alt)) || any(alt != round(alt)) || any(alt < 0))
    stop("Argument `alt` is invalid")
  nInd <- nrow(ref)
  nSnps <- ncol(ref)
  if(nInd != nrow(alt) || nSnps != ncol(alt))
    stop("The dimension of the matrix for the reference reads is not the same the dimension of the matrix for alternate reads")
  if(!is.matrix(parhap) || nrow(parhap) != 4 || ncol(parhap) != nSnps || any(is.na(parhap)) || length(unique(as.vector(parhap))) != 2)
    stop("Argument `parhap` is invalid")
  alleles <- unique(as.vector(parhap))
  major = alleles[1]
  minor = alleles[2]
  OPGP = parHapToOPGP(parhap, major=major, minor=minor)
  if(!is.numeric(iter) || length(iter) != 1 || iter != round(iter) || iter < 1 || !is.finite(iter))
    stop("Argument `iter` is invalid")
  if(!is.numeric(burnin) || length(burnin) != 1 || burnin != round(burnin) || burnin < 1 || !is.finite(burnin))
    stop("Argument `burnin` is invalid")
  if(!is.numeric(chains) || length(chains) != 1 || chains != round(chains) || chains < 1 || !is.finite(chains))
    stop("Argument `iter` is invalid")
  
  if(is.null(initval)){
    set.seed(seed)
    startVal <- replicate(chains, c(rnorm(nSnps-1,-5,1), rnorm(nSnps,-5,0.8),
                                    rnorm(1,-5,1),rnorm(1,-5,0.8),rlnorm(1, sdlog=0.5),rlnorm(1, sdlog=0.5)), simplify=F)
  }
  else if(!is.list(initval) || length(initval) != chains || all(lapply(initval, length) == (2*nSnps+3)))
    stop("Argument `initval` is invalid")
  else startVal = initval
  ## run the MH algorithm
  doParallel::registerDoParallel(cores=cores)
  out <- foreach::foreach(chain = 1:chains) %dopar% {
    MH_Bayes_Hir_seq(ref, alt, OPGP, nInd, nSnps, startVal[[chain]], c(burnin,iter), seed*2+13*chain)
  }
  out <- lapply(out, function(x) {
    y = x
    colnames(y) = c(paste0("\u03C1",1:(nSnps-1)),paste0("\u03B5",1:nSnps),
                    paste0("\u03BC",1:2),paste0("\u03C3",1:2))
    return(y)
  })
  names(out) <- paste0("chain",1:length(out))
  return(out)
}