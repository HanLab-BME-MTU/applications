function [x,y]=geomDeformedRect(bs,s)
global globalCurveL globalCurveT globalCurveR globalCurveB;

% Run along the whole curves, and determine the curve length up to each
% point and calculate the arc length for each point:

[arcLengthL , numPtsL]=curveLength(globalCurveL);
[arcLengthT , numPtsT]=curveLength(globalCurveT);
[arcLengthR , numPtsR]=curveLength(globalCurveR);
[arcLengthB , numPtsB]=curveLength(globalCurveB);

% Create a single curve for the whole boundary.
% It is necessary to skip the corners that appear twice, despite the end
% point!
totalCurve=vertcat(globalCurveL(1:end-1,:), globalCurveT(1:end-1,:), globalCurveR(1:end-1,:), globalCurveB);

% Do the same for the whole boundary:
[arcLengthTotal,numPtsTotal]=curveLength(totalCurve);

% The arc length of each point is then given by:
b=2*pi*arcLengthTotal/arcLengthTotal(end);

%**************************************************************************
% In the following, the different results are calculated. There are three 
% cases to be treated: 
% 1) nargin=0; 
% 2) nargin=1;
% 3) nargin=2;
%**************************************************************************

% Per definition, the number of boundaries is 4.
% This value is retured when 1) nargin=0.
nbs=4; 
 
if nargin==0  
    x=nbs;   
    return; 
end

% In case 2) nargin=1, information about the arc-length intervall of each
% edge has to be returned:

% The arc-length intervall of each edge is given by:
startCurveL=2*pi*arcLengthL(1)  /arcLengthTotal(end); %should be 0
endCurveL  =2*pi*arcLengthL(end)/arcLengthTotal(end);

startCurveT=endCurveL+2*pi*arcLengthT(1)  /arcLengthTotal(end);
endCurveT  =endCurveL+2*pi*arcLengthT(end)/arcLengthTotal(end);

startCurveR=endCurveT+2*pi*arcLengthR(1)  /arcLengthTotal(end);
endCurveR  =endCurveT+2*pi*arcLengthR(end)/arcLengthTotal(end);

startCurveB=endCurveR+2*pi*arcLengthB(1)  /arcLengthTotal(end);
endCurveB  =endCurveR+2*pi*arcLengthB(end)/arcLengthTotal(end); %should be 2*pi

dl=[  startCurveL   startCurveT    startCurveR   startCurveB
      endCurveL     endCurveT      endCurveR     endCurveB
            1             1               1               1
            0             0               0               0];

if nargin==1   
    x=dl(:,bs);   
    return;
end

% This following check is necessary since the comsol algorithm meshinit sometimes
% calls the geometry M-file with values s<0 or s>2pi. Although this 
% overflow in s is marginal it causes NaN
% values when the interpolation is done and later on causes meshinit to crash.
% The algorithm initmesh from Matlab doesn't have this problem.

negIndex=(s<0);
while sum(sum(negIndex))>0
     s=s+negIndex*2*pi;
     negIndex=(s<0);
     display('corrected negative values')
end

exIndex=(s>2*pi);
while sum(sum(exIndex))>0
     s=s-exIndex*2*pi;
     exIndex=(s<0);
     display('corrected exceeding values')
end

s_new=zeros(size(s));
x=zeros(size(s));
y=zeros(size(s));

[m,~]=size(s);
for i=1:m
    s_new(i,:)=pdearcl(b',totalCurve',s(i,:),0,2*pi);
    % In principle one could use s instead of s_new here and skip the
    % previous line. This might be helpful for more general problems.
    x(i,:) = interp1(b,totalCurve(:,1),s_new(i,:),'linear');
    y(i,:) = interp1(b,totalCurve(:,2),s_new(i,:),'linear');
end

function [arcLength,numPts]=curveLength(curve)
    numPts=length(curve);
    arcLength=zeros(numPts,1);
    for i=1:numPts
        if i>1
            arcLength(i)=arcLength(i-1)+sqrt((curve(i,1)-curve(i-1,1))^2+(curve(i,2)-curve(i-1,2))^2);
        end
    end
end

end