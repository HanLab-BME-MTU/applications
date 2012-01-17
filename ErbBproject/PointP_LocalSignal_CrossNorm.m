%Calculates local signal in point point process in a second point process
%Where the points in posB is normalized by the number of posA point that
%will count them.
%
%
%Jeffrey L. Werbin
%Harvard Medical School
%
%Last Update: 1/13/2012
%
%
%

function [out]=PointP_LocalSignal_CrossNorm(posA,posB,rA,rB)
%PointP_LocalSignal_CrossNorm takes two related point processes and evaulates
%two measures around each point in posA: the number of posB less than rB
%from posA(i) and the number of other posA less than rA from posA(i)
%both processes are assumed to fill the retangular area defined by the
%maxium and minimum points in the combined position lists. 
%This version of the code normalizes 


%edge = 1-(4/(3*pi))*((r/L)+(r/W))+((11/(3*pi))-1)*(r^2/(L*W));
edge=1;

%The number of total points
n = size(posA,1);
m = size(posB,1);

pos = [posA ; posB]; %concatinating the two lists allows for processing with one dist matrix

%find boundaries
limX = [max(pos(:,1)), min(pos(:,1))];
limY = [max(pos(:,2)), min(pos(:,2))];
W = limX(1)-limX(2);
L = limY(1)-limY(2);

%to analyze only points > max(rA,rB) from the boundries of the rectangle
%removes edge effect problems

Rmax = max(rA,rB);
x = posA(:,1)-limX(2);
y = posA(:,2)-limY(2);

ind = find( (x > Rmax) & (x < W -Rmax) & (y > Rmax) & (y < L-Rmax));

%Make distance array
dis = squareform(pdist(pos));

%Uncomputed values are returned as NaNs
out = NaN(n,2);

%computes the normalization factors for all activated points
norm = dis(n+1:end,1:n)<= rB & dis(n+1:end,1:n) >0;
norm = 1./sum(norm,2); %yielding the number of points in posA less than rB from each posB signal

%sanitizes norm of infinities
norm(isinf(norm))=0;

for i=ind'
    %selects one point in a
    tempA = dis(i,1:n);
    tempB = dis(i,n+1:end);
    
    %applies the local conditions note that it specifically removes the
    %zero point to exclude it from the count
    tempA = tempA<rA & tempA>0;
    tempB = tempB<rB & tempB>0;
    
    %Normalization, vector math takes the sum
    tempB = tempB*norm;
    
    out(i,1)= tempB/(edge*pi*rB^2);
    out(i,2)= sum(tempA)/(edge*pi*rA^2);
end

end
    