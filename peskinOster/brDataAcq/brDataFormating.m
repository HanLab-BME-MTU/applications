function [mt1Vel,cAngle1T,mt2Vel,cAngle2T]=brDataFormating(idlisttrack,dataProperties,lastResult,stateWanted);
%BRFORMATING put the data in the needed format for the optimzation

%SYNOPSIS
%[mt1Vel,cAngle1T,mt2Vel,cAngle2T]=brDataFormating(idlisttrack,dataProperties,lastResult,stateWanted)
%    
%
%INPUT 
%       idlisttrack     : list from jonas tracking program
%       dataProperties  : properties from jonas tracking program
%       lastResult      : last result from jonas tracking program
%       sateWanted      : sting, GG growth-growth case, SS SG or GS
%OUTPUT
%       mt1Vel          : vector with the velocity for first MT
%       mt2Vel          : vector with the velocity for second MT
%       cAngle1T        : vector with the cosinus of the angle between the
%                         first MT and the inter-centromere axis
%       cAngle2T        : vector with the cosinus of the angle between the
%                         second MT and the inter-centromere axis
%COMMENTS: Not anymore used in the latest version



[mt1State,mt2State,mt1Velocity,mt2Velocity,cAngle1T,cAngle1TdT,cAngle2T,cAngle2TdT]=brGetAnaData(idlisttrack,dataProperties,lastResult);

switch stateWanted
    case 'GG'
        index=find(mt1State>0&mt2State>0);
    case 'SS'
        index=find(mt1State<0&mt2State<0);
    case 'SG'
        index=find(mt1State<0&mt2State>0);
    case 'GS'
        index=find(mt1State>0&mt2State<0);
end

mt1Vel=mt1Velocity(index);
cAngle1T=cAngle1T(index);

mt2Vel=mt2Velocity(index);
cAngle2T=cAngle2T(index)

if length(mt1Vel)~=length(mt2Vel)
    disp('dd');
end
if length(mt1Vel)~=length(cAngle2T);
    disp('ff');
end


