## Function for converting parental haplotypes to OPGPs
parHapToOPGP <- function(parHap, major="A", minor="B"){
  return (apply(parHap,2,function(x){
    if (all(x==major)){
      return(13) 
    } else if (x[1]==x[2] & x[3]==x[4] & x[1]==major & x[3]==minor){
      return(14)
    } else if (x[1]==x[2] & x[3]==x[4] & x[3]==major & x[1]==minor){
      return(15)
    } else if (all(x==minor)){
      return(16)
    } else if(sum(x==major)==2){
      return((x[1]==minor)+(x[3]==minor)*2+1)
    } else if (all(x[3:4]==major)){
      return((x[1]==minor)+5)
    } else if (all(x[1:2]==major)){
      return((x[3]==minor)+9)
    } else if (all(x[3:4]==minor)){
      return((x[1]==minor)+7)
    } else if (all(x[1:2]==minor)){
      return((x[3]==minor)+11)
    } 
  }) 
  )}

