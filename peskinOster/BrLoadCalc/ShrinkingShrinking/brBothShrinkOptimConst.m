function [c,ceq]=brBothShrinkOptimConst(x,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);
%brBothShrinkOptimConst return the const. to optimize in both shrink MT case 

switch pFixe

    case 'free'
		param1=x(1:3);
		param2=x(4:6);
		
		% discretization
		domega=0.001;
		omegaDiscret=[0:domega:10];
		velDiscret1=brShrinkDirect(param1,omegaDiscret,delta);
		velDiscret2=brShrinkDirect(param2,omegaDiscret,delta);
		velMax1=max(velDiscret1);
		velMin1=min(velDiscret1);
		velMax2=max(velDiscret2);
		velMin2=min(velDiscret2);
		
		velData1=abs(velData1);
		velData2=abs(velData2);
		
		c(1)=0.8*velMin1-min(velData1);
		c(2)=0.8*velMin2-min(velData2);
		
		c(3)=0.8*max(velData1)-velMax1;
		c(4)=0.8*max(velData2)-velMax2;
		
		c(5)=0.8*velMax2-max(velData2);
		c(6)=0.8*velMax1-max(velData1);
		c(7)=-x(1);
		c(8)=-x(2);
		c(9)=-x(3);
		c(10)=-x(4);
		c(11)=-x(5);
		c(12)=-x(6);
        c=[c -velDiscret1 -velDiscret2];

		
		% c(7)=(1/(2*x(3))*(1+((1+x(3))/(1-x(3)))^0.5))-x(2)/x(1);
		% c(8)=(1/(2*x(6))*(1+((1+x(6))/(1-x(6)))^0.5))-x(5)/x(4);
		
		
		ceq=[];
        
    case 'fixed'
        
                
        p1=(x(1)-((delta*x(2)*x(1))/vUnload))/x(2)+1;
        p2=(x(3)-((delta*x(4)*x(3))/vUnload))/x(4)+1;
        
  		param1=[x(1:2) p1];
		param2=[x(3:4) p2];
		
		% discretization
		domega=0.001;
		omegaDiscret=[0:domega:10];
		velDiscret1=brShrinkDirect(param1,omegaDiscret,delta);
		velDiscret2=brShrinkDirect(param2,omegaDiscret,delta);
		velMax1=max(velDiscret1);
		velMin1=min(velDiscret1);
		velMax2=max(velDiscret2);
		velMin2=min(velDiscret2);
		
		velData1=abs(velData1);
		velData2=abs(velData2);
		
		c(1)=0.8*velMin1-min(velData1);
		c(2)=0.8*velMin2-min(velData2);
		
		c(3)=0.80*max(velData1)-velMax1;
		c(4)=0.80*max(velData2)-velMax2;
		
		c(5)=0.80*velMax2-max(velData2);
		c(6)=0.80*velMax1-max(velData1);
		c(7)=-x(1);
		c(8)=-x(2);
		c(9)=-x(3);
		c(10)=-x(4);
        c=[c -velDiscret1 -velDiscret2];
		

		
		% c(7)=(1/(2*x(3))*(1+((1+x(3))/(1-x(3)))^0.5))-x(2)/x(1);
		% c(8)=(1/(2*x(6))*(1+((1+x(6))/(1-x(6)))^0.5))-x(5)/x(4);
		
		
		ceq=[];      
end