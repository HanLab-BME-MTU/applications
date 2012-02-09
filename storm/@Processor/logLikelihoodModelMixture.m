function logL = logLikelihoodModelMixture(logL1,logL2,n1,n2)

n = n1+n2;
logL = n1*log(n1/n)+logL1+n2*log(n2/n)+logL2;

end

