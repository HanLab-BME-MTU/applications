function [c,ceq]=brDirectConst(x,velData1,velData2,cAngle1,cAngle2,delta,v0,caseOptim,multi,pFixe)




% inptut check

if (length(velData1)~=length(velData2) | length(velData1)~=length(cAngle1) | length(velData1)~=length(cAngle2))
    error('All the data input have not the same length');
end

if delta~=8e-3 & ~strcmp(caseOptim,'GG')
    warning(' Your delta is not 8e^{-3} \mum');
end

% put all the vector in column form

if size(velData1)~=[length(velData1) 1]
    velData1=velData1';
    velData2=velData2';
end
if size(cAngle1)~=[length(cAngle1) 1]
    cAngle1=cAnlge1';
    cAngle2=cAngle2';
end

if size(x)~=[length(x) 1]
    x=x';
end


switch pFixe
    case 'free'
		
		if strcmp(caseOptim,'GG') & (length(x)-length(velData1))==4
            
			paramG1=x(1:2);
			paramS1=[NaN NaN NaN];
			paramG2=x(3:4);
			paramS2=[NaN NaN NaN];
            
            omegaEffectif  =x(5:end);
		elseif strcmp(caseOptim,'SS') & (length(x)-length(velData1))==6
			paramG1=[NaN NaN];
			paramS1=x(1:3);
			paramG2=[NaN NaN];
			paramS2=x(4:6);
            
            omegaEffectif  =x(7:end);    
            
		elseif strcmp(caseOptim,'GS') & (length(x)-length(velData1))==5
			paramG1=x(1:2);
			paramS1=[NaN NaN NaN];
			paramG2=[NaN NaN];
			paramS2=x(3:5);
            
            omegaEffectif  =x(6:end);  
		
        elseif strcmp(caseOptim,'SG') & (length(x)-length(velData1))==5
			paramG1=[NaN NaN];
			paramS1=x(1:3);
			paramG2=x(4:5);
			paramS2=[NaN NaN NaN];
            
            omegaEffectif  =x(6:end); 
            
        elseif strcmp(caseOptim,'All') & (length(x)-length(velData1))==10
            paramG1=x(1:2);
            paramS1=x(3:5);
            paramG2=x(6:7);
            paramS2=x(8:10);
            
            omegaEffectif=x(11:end);
		else
            error('The input x correspond not with the choosen case');
		end
        
    case 'fixed'
        if isempty(v0)
            error('You choose the free option, but you enter [] for the unload velocity');
        end
        
		if strcmp(caseOptim,'GG') & (length(x)-length(velData1))==4
            
			paramG1=x(1:2);
			paramS1=[NaN NaN NaN];
			paramG2=x(3:4);
			paramS2=[NaN NaN NaN];
            
            omegaEffectif  =x(5:end);
		elseif strcmp(caseOptim,'SS') & (length(x)-length(velData1))==4
            p1=(x(1)-((delta*x(2)*x(1))/v0))/x(2)+1;
            p2=(x(3)-((delta*x(4)*x(3))/v0))/x(4)+1;
			paramG1=[NaN NaN];
			paramS1=[x(1:2); p1];
			paramG2=[NaN NaN];
			paramS2=[x(3:4) ;p2];
            
            omegaEffectif  =x(5:end);    
            
		elseif strcmp(caseOptim,'GS') & (length(x)-length(velData1))==4
            
            p2=(x(3)-((delta*x(4)*x(3))/v0))/x(4)+1;
			paramG1=x(1:2);
			paramS1=[NaN NaN NaN];
			paramG2=[NaN NaN];
			paramS2=[x(3:4); p2];
            
            omegaEffectif  =x(5:end);  
		
            elseif strcmp(caseOptim,'SG') & (length(x)-length(velData1))==4
                
            p1=(x(1)-((delta*x(2)*x(1))/v0))/x(2)+1;
			paramG1=[NaN NaN];
			paramS1=[x(1:2); p1];
			paramG2=x(3:4);
			paramS2=[NaN NaN NaN];
            
            omegaEffectif  =x(5:end); 
            
        elseif strcmp(caseOptim,'All') & (length(x)-length(velData1))==8
            p1=(x(3)-((delta*x(4)*x(3))/v0))/x(4)+1;
            p2=(x(7)-((delta*x(8)*x(7))/v0))/x(8)+1;
            paramG1=x(1:2);
 			paramS1=[x(3:4); p1];
			paramG2=x(5:6);
			paramS2=[x(7:8); p2]; 
            
            omegaEffectif=x(9:end);
            
		else
            error('The input x correspond not with the choosen case');
		end
end
        
        

% calculation of the real omegas
% from the cos(theta) and the effective omega

omega1 = omegaEffectif./cAngle1;
omega2 = omegaEffectif./cAngle2;

% determiation of the mt state

indexMt1Growth=find(velData1>0);
indexMt1Shrink=find(velData1<0);

indexMt2Growth=find(velData2>0);
indexMt2Shrink=find(velData2<0);

velData1=abs(velData1);
velData2=abs(velData2);





% introduction of  state variable commun for both
% mts: 1 if mt1 G 2 if mt1 S and 10 if mt2 G 20 if mt2 S
% global variable: GG 11 / SS 22 / GS 21 / SG 12

stateMt1Tmp(indexMt1Growth)=1;
stateMt1Tmp(indexMt1Shrink)=2;
stateMt2Tmp(indexMt2Growth)=10;
stateMt2Tmp(indexMt2Shrink)=20;

state=stateMt1Tmp + stateMt2Tmp;

%  G

v1=peskin(paramG1,omega1(indexMt1Growth),delta);
v2=peskin(paramG2,omega2(indexMt2Growth),delta);
if ~isempty(v1) |~isempty(v2)
	omegaEff1=x(end-length(velData1)+1:end)./cAngle1;
	omegaEff2=x(end-length(velData1)+1:end)./cAngle2;
    deltaOmega11=omega1-1.2*peskinInv(velData1,x(1:2),delta);
    deltaOmega12=omega2-1.2*peskinInv(velData2,x(3:4),delta);
    deltaOmega21=(peskinInv(velData1,x(1:2),delta)-1.2*omega1);
    deltaOmega22=(peskinInv(velData2,x(3:4),delta)-1.2*omega2);
	cG= [-v1; -v2; -omegaEff1 ; omegaEff1-5;-omegaEff2 ; omegaEff2-5;deltaOmega11;deltaOmega12;deltaOmega21;deltaOmega22];
	cG= [cG ;(paramG1(1)*ones(size(omega1))-max(v1)*1.15./(delta*(1-exp(-omega1)))); (paramG2(1)*ones(size(omega2))-max(v2)*1.15./(delta*(1-exp(-omega2))))]; 
	clear v1 v2;
else
    cG=[];
end

 %  S
 
v1=brShrinkDirect(paramS1,omega1(indexMt1Shrink),delta);
v2=brShrinkDirect(paramS2,omega2(indexMt2Shrink),delta);

cS=[-v1 ;-v2];
if ~isempty(v1)
    cSTmp(1)=0.80*max(velData1)-max(v1);
else 
    cSTmp(1)=-1;
end
if ~isempty(v2)
	cSTmp(2)=0.80*max(velData2)-max(v2);
else
    cSTmp(2)=-1;
end
    
cS=[cS; cSTmp'];

c=[cG ;cS];

ceq=[];

