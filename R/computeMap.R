##########################################################################
# Bayesian genotyping Uncertainty with Sequencing data and linkage MAPping (BUSMap)
# Copyright 2019 Timothy P. Bilton <tbilton@maths.otago.ac.nz>
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
#' Construct linkage map using a Bayesian hierarchical hidden Markov model (HMM)
#' 
#' Function that implements the Metropolis-Hastings (MH) algorithm for obtaining posterior samples 
#' of the Bayesian hierarchical hidden Markov model (HMM) for linkage maps in full-sib families 
#' with high-throughput sequencing data.
#' 
#' @usage
#' computeMap(ref, alt, OPGP, initval=NULL, iter=30000, burnin=5000, chains=3, seed=1, cores=chains)
#' 
#' @param ref Non-negative integer matrix containing the read counts for the reference allele
#' @param alt Non-negative integer matrix containing the read counts for the alternate allele
#' @param OPGP Integer vector specifying the OPGPs 
#' @param initval List giving the starting values for the MH algorithm value
#' @param iter Integer value of the number of iterations in the MH algorithm (excluding burn-in period) for each chain
#' @param burnin Integer value of the number of interations in the burn-in/adaptive phase of the MH algorithm for each chain
#' @param chains Integer value giving the number of parallel chains to use  
#' @param seed Numeric value giving the seed.
#' @param cores Integer value giving the number of cores for parallelization of the MH chains.
#' Must be non-negative and no more than \code{chains}. 
#' 
#' @name computeMap
#' @author Timothy P. Bilton
#' @examples
#' 
#' data(manuka)
#' ref <- manuka$ref
#' alt <- manuka$alt
#' OPGP <- manuka$OPGP
#' 
#' computeMap(ref, alt, OPGP)
#'
#' @export computeMap
computeMap <- function(ref, alt, OPGP, initval=NULL, iter=30000, burnin=5000, chains=3, seed=1, cores=chains){
  ## Do some check
  if(!is.matrix(ref) || !is.numeric(ref) || any(is.na(ref)) || any(ref != round(ref)) || any(ref < 0))
    stop("Argument `ref` is invalid")
  if(!is.matrix(alt) || !is.numeric(alt) || any(is.na(alt)) || any(alt != round(alt)) || any(alt < 0))
    stop("Argument `alt` is invalid")
  nInd <- nrow(ref)
  nSnps <- ncol(ref)
  if(nInd != nrow(alt) || nSnps != ncol(alt))
    stop("The dimension of the matrix for the reference reads is not the same the dimension of the matrix for alternate reads")
  if(!is.vector(OPGP) || any(is.na(OPGP)) || !is.integer(OPGP) || OPGP < 1 || OPGP > 12)
    stop("Argument `OPGP` is invalid")
  if(!is.numeric(iter) || length(iter) != 1 || iter != round(iter) || iter < 1 || !is.finite(iter))
    stop("Argument `iter` is invalid")
  if(!is.numeric(burnin) || length(burnin) != 1 || burnin != round(burnin) || burnin < 1 || !is.finite(burnin))
    stop("Argument `burnin` is invalid")
  if(!is.numeric(chains) || length(chains) != 1 || chains != round(chains) || chains < 1 || !is.finite(chains))
    stop("Argument `iter` is invalid")
  
  if(is.null(initval))
    startVal <- replicate(chains, c(rnorm(nSnps-1,-5,1), rnorm(nSnps,-5,0.8),
                                    rnorm(1,-5,1),rnorm(1,-5,0.8),rlnorm(1, sdlog=0.5),rlnorm(1, sdlog=0.5)), simplify=F)
  else if(!is.list(initval) || length(initval) != chains || all(lapply(initval, length) == (2*nSnps+3)))
    stop("Argument `initval` is invalid")
  else startVal = initval
  ## run the MH algorithm
  doParallel::registerDoParallel(chains)
  out <- foreach::foreach(chain = 1:chains) %dopar% {
    MH_Bayes_Hir(ref, alt, OPGP, nInd, nSnps, startVal[[chain]], c(burnin,iter), seed+13*chain)
  }
  out <- lapply(out, function(x) {
    y = x
    colnames(y) = c(paste0("\u03C1",1:(nSnps-1)),paste0("\u03B5",1:nSnps),
                    paste0("\u03BC",1:2),paste0("\u03C3",1:2))
    return(y)
  })
  return(out)
}