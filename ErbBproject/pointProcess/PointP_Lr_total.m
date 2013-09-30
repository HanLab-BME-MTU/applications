function [kr,lr]=PointP_Lr_total(pos,r)
%PointP_Lr calculates Besags L(r) statistic with Ward & Ferandino's edge
%correction method (1999) that is global and thus fast.
%Ward, J. S. and Ferrandino, F. J. 1999. New derivation reduces bias and
%increases power of Ripley’s L index. / Ecol. Modell. 116: 225/236.
%r is a range
%
%Calculates L(r) second order statistic for a range of r for all
%points 
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 9/6/2011
tic
W = max(pos(:,1))-min(pos(:,1));
L = max(pos(:,2))-min(pos(:,2));
A= W*L;

%edge=1;

dis = squareform(pdist(pos));


n = size(pos,1);
k = size(r,2);
kr = zeros(k,1);
lr = zeros(k,1);

for i=1:k
    temp = dis<r(i) & dis>0;
    %Ward Ferrandino correction for rectangular sampling domain
    edge = 1-(4/(3*pi))*((r(i)/L)+(r(i)/W))+((11/(3*pi))-1)*(r(i)^2/(L*W));
    kr(i) = sum(sum(temp))*((A)/((n-1)*n*edge));
    lr(i)= sqrt(kr(i)/pi)-r(i);
end
toc
end
    