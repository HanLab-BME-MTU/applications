%Calculates L(r) second order statistic for a particular value of r for all
%points 
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 9/6/2011
%
%
%

function [lr]=PointP_Lr(pos,r)
%PointP_Lr calculates Besags L(r) statistic with Ward & Ferandino's edge
%correction method (1999) that is global and thus fast.
%Ward, J. S. and Ferrandino, F. J. 1999. New derivation reduces bias and
%increases power of Ripley’s L index. / Ecol. Modell. 116: 225/236.

W = max(pos(:,1))-min(pos(:,1));
L = max(pos(:,2))-min(pos(:,2));

%edge = 1-(4/(3*pi))*((r/L)+(r/W))+((11/(3*pi))-1)*(r^2/(L*W));
edge=1;

dis = squareform(pdist(pos));
dis = dis<r & dis>0;

n = size(pos,1);
lr = zeros(n,1);

for i=1:n
    lr(i)= sqrt(sum(dis(i,:))/(edge*pi))-r;
end

end
    