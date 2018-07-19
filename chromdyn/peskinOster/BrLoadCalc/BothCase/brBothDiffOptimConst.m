function [c,ceq]=brBothDiffOptimConst(x,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);
% brBothDiffOptimConst generate the constraint for the growth-shrink or
% shrink-grwth case
% INPUT
%     x         :  paramters of pol/depol, input in the correct way: GS first pol
%                  then depol for SG first depol then pol
%     velData1  :  velocity of the MT1 with sign if GS>0 if SG <0
%     velData2  :  velocity of MT2 with sign
%     cAngle1   :  cosines of angle between MT1 and inter-centromere axis
%     cAngle2   :  cosines of angle between MT2 and inter-centromere axis
%     delta     :  1/13 of subunit length
%     vUnload   :  unload velocity of the shrink case
%     pFixe     :  fixed or free regarding if you optimize with p or
%     without
% OUTPUT
%     c         : conatraint inequality in the form for fminimax or fmincon 

% COMMENT : the free case is not always update
if length(find(velData1<0))==length(velData1) & length(find(velData2>0))==length(velData2)
    state='SG';
elseif length(find(velData1>0))==length(velData1) & length(find(velData2<0))==length(velData2)
    state='GS';
else
    error('You do not input velocity of growing/shrinkage or shrionkage/growing case');
end


domega=0.001;
omegaDiscret=[0:domega:10];

switch pFixe
    case 'free'
		
		switch state
            case 'SG'
                param1=x(1:3);
                
                % discretization
		
				velDiscret1=brShrinkDirect(param1,omegaDiscret,delta);
				
				velMax1=max(velDiscret1);
				velMin1=min(velDiscret1);
				
				velData1=abs(velData1);
                
                c(1)=0.6*velMin1-min(velData1);
				c(2)=0.6*max(velData1)-velMax1;
				c(3)=0.6*velMax1-max(velData1);
                
                alpha2=x(4);
				beta2=x(5);
			
				omega=[0:0.1:5];
				v2=delta/13*(alpha2*exp(-omega)-beta2);
                v2Max=delta/13*(alpha2-beta2);
                
                c=[c -v2 ];
				
		
		
				
                
            case 'GS'
                
                param2=x(3:5);
                
                velDiscret2=brShrinkDirect(param2,omegaDiscret,delta);
                velMax2=max(velDiscret2);
				velMin2=min(velDiscret2);
                
                velData2=abs(velData2);
                
                c(1)=0.6*velMin2-min(velData2);
                c(2)=0.6*max(velData2)-velMax2;
                c(3)=0.6*velMax2-max(velData2);
                
                alpha1=x(1);
				beta1=x(2);
			
				omega=[0:0.1:5];
                v1=delta/13*(alpha1*exp(-omega)-beta1);
                v1Max=delta/13*(alpha1-beta1);
                
                c=[c -v1];
                
		end
		
		ceq=[];
                
            case 'fixed'
		
		switch state
            case 'SG'
                
                p1=(x(1)-((delta*x(2)*x(1))/vUnload))/x(2)+1;
                param1=[x(1:2) p1];
                
                % discretization
		
				velDiscret1=brShrinkDirect(param1,omegaDiscret,delta);
				
				velMax1=max(velDiscret1);
				velMin1=min(velDiscret1);
				
				velData1=abs(velData1);
                
                c(1)=velMin1-min(velData1);
				c(2)=max(velData1)-velMax1;
				c(3)=0.6*velMax1-max(velData1);
                
                alpha2=x(3);
				beta2=x(4);
			
				omega=[0:0.1:5];
				v2=delta/13*(alpha2*exp(-omega)-beta2);
                v2Max=delta/13*(alpha2-beta2);
                c(4)=0.6*v2Max-max(velData2);
                
                c=[c -v2];
				
		
		
				
                
            case 'GS'
                
                p2=(x(3)-((delta*x(4)*x(3))/vUnload))/x(4)+1;
                
                param2=[x(3:4) p2];
                
                velDiscret2=brShrinkDirect(param2,omegaDiscret,delta);
                velMax2=max(velDiscret2);
				velMin2=min(velDiscret2);
                
                velData2=abs(velData2);
                
                c(1)=velMin2-min(velData2);
                c(2)=max(velData2)-velMax2;
                c(3)=0.6*velMax2-max(velData2);
                
                alpha1=x(1);
				beta1=x(2);
			
				omega=[0:0.1:5];
                v1=delta/13*(alpha1*exp(-omega)-beta1);
                v1Max=delta/13*(alpha1-beta1);
                
                c(4)=0.6*v1Max-max(velData1);
                
                c=[c -v1];
                
		end
		
		ceq=[];
end