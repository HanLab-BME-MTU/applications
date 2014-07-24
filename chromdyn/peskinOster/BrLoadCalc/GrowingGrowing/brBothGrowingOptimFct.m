function [fct]=brBothGrowingOptimFct(x,velData1,velData2,cAngle1,cAngle2,delta);
% brBothGrowingOptimFct generate the obj. fct for the growth-growth case
% INPUT
%     x         :  paramters of pol
%     velData1  :  velocity of the MT1 w>0
%     velData2  :  velocity of MT2 >0
%     delta     :  1/13 of subunit length
%     without
% OUPUT 
%     fct       :  objective function
% COMMENT : the free case is not always update


        

alpha1=x(1);
beta1=x(2);
alpha2=x(3);
beta2=x(4);

% calculation of omega
omega1=peskinInv(velData1,x(1:2),delta);
omega2=peskinInv(velData2,x(3:4),delta);



% fct= (delta(proj(omegaVect on centtomer Line)))^2

if size(cAngle1)~=size(omega1)
    cAngle1=cAngle1';
end
if size(cAngle2)~=size(omega2)
    cAngle2=cAngle2';
end

fct=(cAngle1.*omega1-cAngle2.*omega2).^2;


		
