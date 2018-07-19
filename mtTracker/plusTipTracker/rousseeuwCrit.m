function [outliers]=rousseeuwCrit(data)
% rousseeuwCrit finds outliers using Rousseeuw's criterion (ref. Gaudenz)
%
% SYNOPSIS: [outliers]=rousseeuwCrit(data)
%
% data    : n x m matrix where outliers are detected separately from the m
%           populations of n observations
% outliers: n x m matrix where outliers are 1 and inliers are 0


m=nanmedian(data);
m=repmat(m,size(data,1),1);
dist2med=abs(data-m);
So=1.4826*sqrt(nanmedian(dist2med.^2));
So=repmat(So,size(data,1),1);
outIdx=find(dist2med>3.*So);
outliers=zeros(size(data));
outliers(outIdx)=1;