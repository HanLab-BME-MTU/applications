function [mt1Vel,cAngle1T,mt2Vel,cAngle2T]=brDataFormating(idlisttrack,dataProperties,lastResult,stateWanted);





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


