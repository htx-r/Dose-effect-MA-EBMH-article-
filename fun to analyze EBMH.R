# function to compute the relative odds ratio as: odds of non-reference dose / odds of refernce dose using metabin
library(meta)
createORreference.fun=function(r,n)
{

  logOR=c(0)
  selogOR=c(NA)

  for(i in 2:c(length(n)))
  {
    calculate=metabin(r[i],n[i],r[1],n[1],sm="OR")
    logOR=c(logOR,calculate$TE)
    selogOR=c(selogOR,calculate$seTE)

  }
  return(cbind(logOR=logOR,selogOR=selogOR))
}
