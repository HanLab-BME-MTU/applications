function [i,j]=indexij(n,m,k)

% provides the matrix index i,j from the linear index k returned by many matlab routines, given the size(n,m) nrows! of the matrix 
 
 if or(k<1,k>n*m)
   % warning('returning indices out of the matrix range! index will not be valid')
 end
     

 tmp1=rem(k,n);
    
    if tmp1==0.
    i=n;
    else
    i=tmp1;
    end
    
 tmp2=floor(k/n);
  
    if rem(k,n)==0.
    j=tmp2;
    else
    j=tmp2+1;
    end
    