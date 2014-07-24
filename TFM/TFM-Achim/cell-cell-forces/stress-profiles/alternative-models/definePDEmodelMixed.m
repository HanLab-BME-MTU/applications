function [b,c,a,f]=definePDEmodelAllDirichlet

% This is not that nice and should be converted into Matrices:
c1='2*((gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio))))+(2*((gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio)))).*(gfun_poisson_ratio)./(1-(gfun_poisson_ratio)))';
c2='0';
c3='(gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio)))';
c4='0';
c5='(gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio)))';
c6='2*((gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio)))).*(gfun_poisson_ratio)./(1-(gfun_poisson_ratio))';
c7='0';
c8='(gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio)))';
c9='0';
c10='2*((gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio))))+(2*((gfun_youngs_mod(x,y))./(2*(1+(gfun_poisson_ratio)))).*(gfun_poisson_ratio)./(1-(gfun_poisson_ratio)))';

c(1 ,1:size(c1 ,2))=c1;
c(2 ,1:size(c2 ,2))=c2;
c(3 ,1:size(c3 ,2))=c3;
c(4 ,1:size(c4 ,2))=c4;
c(5 ,1:size(c5 ,2))=c5;
c(6 ,1:size(c6 ,2))=c6;
c(7 ,1:size(c7 ,2))=c7;
c(8 ,1:size(c8 ,2))=c8;
c(9 ,1:size(c9 ,2))=c9;
c(10,1:size(c10,2))=c10;


a=['0'
   '0'
   '0'
   '0'];

f =['gfun_vol_f_xcomp(x,y)'                                                                                                                                              
    'gfun_vol_f_ycomp(x,y)'];


b =[ 2     2     2     2
     2     0     0     0
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1    48    48    48
     1    48    48    48
     1    48    48    48
     1    48    48    48
     1    48    48    48
     1    48    48    48
    48    48    48    48
    48    48    48    48
    48    48    48    48
    48    48    48    48
    48    48    48    48
    48    48    48    48
    49    49    49    49
    48    48    48    48
    48    48    48    48
    49    49    49    49
    48    48    48    48
    48    48    48    48];