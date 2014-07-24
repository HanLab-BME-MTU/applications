function [kLims, hLims, names,colors] = getCurvCategories

%Units are in 1/microns and the boundaries are COMPLETELY arbitrary
%(originally defined as STD of given curvature for one cell)

%hSmall = .2;
kSmall = .1;


kLims = [ kSmall Inf;
          kSmall Inf;
         -kSmall kSmall;
         -Inf   -kSmall];
hLims = [0      Inf;
        -Inf    0;
        -Inf    Inf;
        -Inf    Inf];
% 
% funs = {@(x)(K > kSmall & H > hSmall),...
%         @(x)(K > -kSmall & H < 0),...
%         @(x)(K > -kSmall & K < kSmall & H > 0),...
%         @(x)(K < -kSmall & H > -hSmall & H < Inf)};           

    
names = {'k1 pos k2 pos',
         'k2 neg k1 neg',
         'k0 pos k2 0',
         'k1 pos k2 neg'};
colors = [1 1 1 ;
          1 0 1 ;
          0 1 0 ;
          1 1 0 ];              