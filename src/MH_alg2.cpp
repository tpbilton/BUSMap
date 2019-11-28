#include <Rcpp.h>
#include <random>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]



double Qentry(int OPGP,double Kaa,double Kab, double Kbb,int elem){
  switch(OPGP){
  case 1:
    if(elem == 1)
      return Kbb;
    else if ((elem == 2)|(elem == 3))  
      return Kab;
    else if (elem == 4)
      return Kaa;
  case 2:
    if(elem == 3)
      return Kbb;
    else if ((elem == 1)|(elem == 4))
      return Kab;
    else if (elem == 2)
      return Kaa;
  case 3:
    if(elem == 2) 
      return Kbb;
    else if ((elem == 1)|(elem == 4))
      return Kab;
    else if (elem == 3)
      return Kaa;
  case 4:
    if(elem == 4) 
      return Kbb;
    else if ((elem == 2)|(elem == 3))
      return Kab;
    else if (elem == 1)
      return Kaa;
  case 5:
    if ((elem == 1)|(elem == 2))
      return Kab;
    else if ((elem == 3)|(elem == 4))
      return Kaa;
  case 6:
    if ((elem == 1)|(elem == 2))
      return Kaa;
    else if ((elem == 3)|(elem == 4))
      return Kab;
  case 7:
    if ((elem == 1)|(elem == 2))
      return Kbb;
    else if ((elem == 3)|(elem == 4))
      return Kab;
  case 8:
    if ((elem == 1)|(elem == 2))
      return Kab;
    else if ((elem == 3)|(elem == 4))
      return Kbb;
  case 9:
    if ((elem == 1)|(elem == 3))
      return Kab;
    else if ((elem == 2)|(elem == 4))
      return Kaa;
  case 10:
    if ((elem == 1)|(elem == 3))
      return Kaa;
    else if ((elem == 2)|(elem == 4))
      return Kab;
  case 11:
    if ((elem == 1)|(elem == 3))
      return Kbb;
    else if ((elem == 2)|(elem == 4))
      return Kab;
  case 12:
    if ((elem == 1)|(elem == 3))
      return Kab;
    else if ((elem == 2)|(elem == 4))
      return Kbb;
  case 13:
    return Kaa;
  case 14:
    return Kab;
  case 15:
    return Kab;
  case 16:
    return Kbb;
  } // end of Switch
  return -1;
}

// Function for returning a specified enetry of the transition matrix for a given recombination fraction value
double Tmat(int s1, int s2, double rval){
  int sSum = s1 + s2*4;
  if((sSum == 0)|(sSum == 5)|(sSum == 10)|(sSum == 15))
    return (1-rval)*(1-rval);
  else if((sSum == 3)|(sSum == 6)|(sSum == 9)|(sSum == 12))
    return rval*rval;
  else
    return (1-rval)*rval;
}

double inv_logit(double x){
  return R::plogis(x, 0.0, 1.0, 1, 0);
}

double inv_cloglog(double x){
  return (1.0 - exp(-1.0*exp(x)))/2.0;
}

double computell(NumericVector para, int nInd, int nSnps, 
                 IntegerMatrix ref, IntegerMatrix alt, IntegerVector OPGP, 
                 NumericMatrix pAB){
  
  int ind, snp;  
  NumericVector rf(nSnps - 1);
  NumericVector ep(nSnps);
  for(snp=0; snp < nSnps-1; snp++){
    rf[snp] = inv_cloglog(para[snp]);
    ep[snp] = inv_logit(para[nSnps - 1 + snp]);
  }
  ep[nSnps-1] = inv_logit(para[2*nSnps - 2]);
  //Rcpp::Rcout << "rf :" << rf[4] << std::endl; 
  //Rcpp::Rcout << "ep :" << ep[1] << std::endl; 
  
  // define the density values for the emission probs
  NumericMatrix pAA(nInd, nSnps);
  NumericMatrix pBB(nInd, nSnps);
  //  #pragma omp parallel for num_threads(nThreads_c) private(index, snp)
  for(ind = 0; ind < nInd; ind++){
    for(snp = 0; snp < nSnps; snp++){
      pAA(ind, snp) = pow(1.0 - ep[snp], ref(ind, snp)) * pow(ep[snp], alt(ind, snp));
      pBB(ind, snp) = pow(1.0 - ep[snp], alt(ind, snp)) * pow(ep[snp], ref(ind, snp));
    }
  }
  
  // Now compute the likelihood
  //#pragma omp parallel for reduction(+:llval) num_threads(nThreads_c) private(sum, s1, alphaDot, alphaTilde, snp, s2, w_new)
  int s1, s2;
  double sum, alphaDot[4], alphaTilde[4], llval, w_new;
  llval = 0;
  for(ind = 0; ind < nInd; ind++){
    // Compute forward probabilities at snp 1
    sum = 0;
    for(s1 = 0; s1 < 4; s1++){
      alphaDot[s1] = 0.25 * Qentry(OPGP[0], pAA(ind, 0), pAB(ind, 0), pBB(ind, 0), s1+1);
      sum = sum + alphaDot[s1];
    }
    // Scale forward probabilities
    for(s1 = 0; s1 < 4; s1++){
      alphaTilde[s1] = alphaDot[s1]/sum;
    }
    
    // add contribution to likelihood
    llval = llval + log(sum);
    
    // iterate over the remaining SNPs
    for(snp = 1; snp < nSnps; snp++){
      // compute the next forward probabilities for snp \ell
      for(s2 = 0; s2 < 4; s2++){
        sum = 0;
        for(s1 = 0; s1 < 4; s1++){
          sum = sum + Tmat(s1, s2, rf[snp-1]) * alphaTilde[s1];
        }
        alphaDot[s2] = Qentry(OPGP[snp], pAA(ind, snp), pAB(ind, snp), pBB(ind, snp), s2+1) * sum;
      }
      // Compute the weight for snp \ell
      w_new = 0;
      for(s2 = 0; s2 < 4; s2++){
        w_new = w_new + alphaDot[s2];
      }
      // Add contribution to the likelihood
      llval = llval + log(w_new);
      // Scale the forward probability vector
      for(s2 = 0; s2 < 4; s2++){
        alphaTilde[s2] = alphaDot[s2]/w_new;
      }
    }
  }
  return llval;
}

// Prior for transformed mean which has beta(1/2,1/2) on logit scale
double pdist_log_logit(double x){
  return 0.5*x - log(1.0+exp(x));
}

// Prior for transformed mean which has truncated beta(1/2,1/2) on[0,0.5] on logit scale
double pdist_log_logit_jeff(double x){
  return 0.5*x - log(1.0+exp(x)) - 0.5*log(2.0+exp(x));
}

// Prior for transformed mean which has truncated beta(1/2,1/2) on[0,0.5] on cloglog scale
double pdist_log_cloglog_jeff(double x){
  return x - exp(x) - 0.5*log(1.0 - exp(-exp(x))) - 0.5*log(1.0 + exp(-exp(x)));
}


// [[Rcpp::export]]
NumericMatrix MH_Bayes_Hir(IntegerMatrix ref, IntegerMatrix alt, IntegerVector OPGP, int nInd, int nSnps, 
                           NumericVector startVal, IntegerVector simPar, unsigned seed){
  
  // Extract simulation paramters
  int nadapt = simPar[0];
  int nmh = simPar[1];
  int npar = 2*nSnps - 1 + 4;
  double df = 3;
  
  // output matrix
  NumericMatrix samples(nadapt + nmh, npar);
  NumericVector pold = clone(startVal);
  NumericVector ptemp(2*nSnps - 1);
  NumericVector jmp(npar, 0.001);      // variance on the proposal distribution
  // Set up parameter vector for likelihood calculations
  int indx;
  for(indx = 0; indx < nSnps - 1; indx++){
    ptemp[indx] = pold[2*nSnps - 1] + pold[indx]*pold[2*nSnps + 1];
    ptemp[indx + nSnps - 1] = pold[2*nSnps] + pold[indx + nSnps-1]*pold[2*nSnps + 2];
  }
  ptemp[2*nSnps - 1 - 1] = pold[2*nSnps] + pold[2*nSnps - 1 - 1]*pold[2*nSnps + 2];
  //Rcpp::Rcout << "ptemp :" << ptemp << std::endl; 
  //Rcpp::Rcout << "jmp :" << jmp << std::endl; 
  //Rcpp::Rcout << "seed :" << seed << std::endl; 
  
  // Set up the random number generator
  //zigg.setSeed(seed*39); // for zigg.norm function
  std::mt19937 eng(seed);
  std::uniform_real_distribution<double> RNGunif(0.0,1.0);
  std::normal_distribution<double> RNGnorm(0.0,1.0);
  
  // Compute the pAB matrix once
  NumericMatrix pAB(nInd,nSnps);
  int ind, snp;
  for(ind = 0; ind < nInd; ind++){
    for(snp = 0; snp < nSnps; snp++){
      pAB(ind,snp) = powl(0.5,ref(ind,snp) + alt(ind, snp));
    }
  }
  
  // Compute the likelihood
  double llold, llnew, llrat;
  //NumericMatrix matForw(), matBack();
  llold = computell(ptemp, nInd, nSnps, ref, alt, OPGP, pAB);
  //Rcpp::Rcout << "llold :" << llold << std::endl; 
  
  // Perform the adaptive phase
  int iter, par;
  double prat, u, pnew;
  for(iter = 0; iter < nadapt; iter++){
    //Rcpp::Rcout << "iter :" << iter << std::endl;
    for(par = 0; par < npar; par++){
      // simulate proposal from N(mu=pold, sd=jmp)
      if(par < (2*nSnps+1)){
        //Rcpp::Rcout << "par :" << par << std::endl;
        pnew = pold[par] + RNGnorm(eng)*jmp[par];
        //// Compute prior ratio
        // rf's
        if(par < (nSnps - 1)){
          prat = exp(0.5*(-pow(pnew,2)+pow(pold[par],2)));
          ptemp[par] = pold[2*nSnps - 1] + pnew*pold[2*nSnps + 1];
        }
        // epsilon
        else if(par < (2*nSnps - 1)){
          prat = exp(0.5*(-pow(pnew,2)+pow(pold[par],2)));
          ptemp[par] = pold[2*nSnps] + pnew*pold[2*nSnps + 2];
        }
        // mu's
        else{ // if(par < (2*nSnps + 1)){
          // Rcpp::Rcout << "ptemp :" << ptemp << std::endl; 
          if(par == (2*nSnps - 1)){
            prat = exp(pdist_log_cloglog_jeff(pnew) - pdist_log_cloglog_jeff(pold[par]));
            for(indx = 0; indx < nSnps - 1; indx++){
              ptemp[indx] = pnew + pold[indx]*pold[2*nSnps + 1];
            }
          } else{
            prat = exp(pdist_log_logit(pnew) - pdist_log_logit(pold[par]));
            for(indx = nSnps - 1; indx < 2*nSnps - 1; indx++){
              ptemp[indx] = pnew + pold[indx]*pold[2*nSnps + 2];
            }
          }
          //Rcpp::Rcout << "ptemp :" << ptemp << std::endl; 
        }
      }
      // sigma's
      else{
        pnew = exp(log(pold[par]) + RNGnorm(eng)*jmp[par]);
        //Rcpp::Rcout << "ptemp :" << pnew << std::endl; 
        if(par == (2*nSnps + 1)){
          for(indx = 0; indx < nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps - 1] + pold[indx]*pnew;
          }
        } else{
          for(indx = nSnps - 1; indx < 2*nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps] + pold[indx]*pnew;
          }
        }
        //pnew = exp(pnew);
        prat = exp((df+1)/2*(log(1+pow(pold[par], 2)/df) - log(1+pow(pnew, 2)/df))) * pnew/pold[par];
      }
      // compute likelihood ratio
      llnew = computell(ptemp, nInd, nSnps, ref, alt, OPGP, pAB);
      llrat = exp(llnew - llold);
      
      // accept with probability mhrat
      u = RNGunif(eng);
      if(llrat*prat > u){
        pold[par] = pnew;
        llold = llnew;
        jmp[par] *= 1.1;
      } else{
        jmp[par] /= 1.1;
        if(par < (nSnps - 1))
          ptemp[par] = pold[2*nSnps - 1] + pold[par]*pold[2*nSnps + 1];
        else if(par < (2*nSnps - 1))
          ptemp[par] = pold[2*nSnps] + pold[par]*pold[2*nSnps + 2];
        else if((par == (2*nSnps-1)) || (par == (2*nSnps + 1))){
          for(indx = 0; indx < nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps - 1] + pold[indx]*pold[2*nSnps + 1];
          }
        } else{
          for(indx = nSnps - 1; indx < 2*nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps] + pold[indx]*pold[2*nSnps + 2];
          }
        }
      }
      // update the samples
      if(par < (nSnps - 1))
        samples(iter, par) = inv_cloglog(ptemp[par]);
      else if(par < (2*nSnps - 1))
        samples(iter, par) = inv_logit(ptemp[par]);
      else
        samples(iter, par) = pold[par];
    }
  }
  
  // Perform the Metrophis-Hastings phase
  for(iter = nadapt; iter < (nadapt + nmh); iter++){
    for(par = 0; par < npar; par++){
      // simulate proposal from N(mu=pold, sd=jmp)
      //Rcpp::Rcout << "par :" << par << std::endl; 
      if(par < (2*nSnps+1)){
        pnew = pold[par] + RNGnorm(eng)*jmp[par];
        //// Compute prior ratio
        // rf's
        if(par < (nSnps - 1)){
          prat = exp(0.5*(-pow(pnew,2)+pow(pold[par],2)));
          ptemp[par] = pold[2*nSnps - 1] + pnew*pold[2*nSnps + 1];
        }
        // epsilon
        else if(par < (2*nSnps - 1)){
          prat = exp(0.5*(-pow(pnew,2)+pow(pold[par],2)));
          ptemp[par] = pold[2*nSnps] + pnew*pold[2*nSnps + 2];
        }
        // mu's
        else{ // if(par < (2*nSnps + 1)){
          if(par == (2*nSnps - 1)){
            prat = exp(pdist_log_cloglog_jeff(pnew) - pdist_log_cloglog_jeff(pold[par]));
            for(indx = 0; indx < nSnps - 1; indx++){
              ptemp[indx] = pnew + pold[indx]*pold[2*nSnps + 1];
            }
          } else{
            prat = exp(pdist_log_logit(pnew) - pdist_log_logit(pold[par]));
            for(indx = nSnps - 1; indx < 2*nSnps - 1; indx++){
              ptemp[indx] = pnew + pold[indx]*pold[2*nSnps + 2];
            }
          }
        }
      }
      // sigma's
      else{
        pnew = exp(log(pold[par]) + RNGnorm(eng)*jmp[par]);
        if(par == 2*nSnps + 1){
          for(indx = 0; indx < nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps - 1] + pold[indx]*pnew;
          }
        } else{
          for(indx = nSnps - 1; indx < 2*nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps] + pold[indx]*pnew;
          }
        }
        //pnew = exp(pnew);
        prat = exp((df+1)/2*(log(1+pow(pold[par], 2)/df) - log(1+pow(pnew, 2)/df))) * pnew/pold[par];
      }
      // compute likelihood ratio
      llnew = computell(ptemp, nInd, nSnps, ref, alt, OPGP, pAB);
      llrat = exp(llnew - llold);
      
      // accept with probability mhrat
      u = RNGunif(eng);
      if(llrat*prat > u){
        pold[par] = pnew;
        llold = llnew;
      } else{
        if(par < (nSnps - 1))
          ptemp[par] = pold[2*nSnps - 1] + pold[par]*pold[2*nSnps + 1];
        else if(par < (2*nSnps - 1))
          ptemp[par] = pold[2*nSnps] + pold[par]*pold[2*nSnps + 2];
        else if((par == (2*nSnps - 1))||(par == (2*nSnps + 1))){
          for(indx = 0; indx < nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps - 1] + pold[indx]*pold[2*nSnps + 1];
          }
        } else{
          for(indx = nSnps - 1; indx < 2*nSnps - 1; indx++){
            ptemp[indx] = pold[2*nSnps] + pold[indx]*pold[2*nSnps + 2];
          }
        }
      }
      // update the samples
      if(par < (nSnps - 1))
        samples(iter,par) = inv_cloglog(ptemp[par]);
      else if(par < (2*nSnps - 1))
        samples(iter, par) = inv_logit(ptemp[par]);
      else
        samples(iter, par) = pold[par];
    }
  }
  return samples;
}

