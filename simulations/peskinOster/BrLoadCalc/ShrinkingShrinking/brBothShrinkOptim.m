function [fct]=brBothShrinkOptim(x,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);
%brBothShrinkOptim return the fct to optimize in both shrink MT case 





switch pFixe
    
    case 'free'
		c=0.1;
		
		
		% discretization
		domega=0.001;
		omegaDiscret=[0:domega:10];
		
		velDiscret1=brShrinkDirect(x(1:3),omegaDiscret,delta);
		velDiscret2=brShrinkDirect(x(4:6),omegaDiscret,delta);
		velMax1=max(velDiscret1);
		velMin1=min(velDiscret1);
		velMax2=max(velDiscret2);
		velMin2=min(velDiscret2);
		
		
		%looking for the shrink velocties in the good range
		
		
		indiceVelOK=find(abs(velData1)>=velMin1 & abs(velData1)<=velMax1 & abs(velData2)>=velMin2 & abs(velData2)<=velMax2 & velData1<0 & velData2<0 ); 
		indiceVelNo=find(abs(velData1)<velMin1  | abs(velData1)>velMax1  | abs(velData2)<velMin2  | abs(velData2)>velMax2 |  velData1>0 | velData2>0 );
		
              
		
		
		velData1=abs(velData1);
		velData2=abs(velData2);
		
		% load calculation 
		
		[omega1(indiceVelOK),omegaHigh1(indiceVelOK),omegaLow1(indiceVelOK)] = brInvShrinkSpline(velData1(indiceVelOK),velDiscret1,omegaDiscret);
		[omega2(indiceVelOK),omegaHigh2(indiceVelOK),omegaLow2(indiceVelOK)] = brInvShrinkSpline(velData2(indiceVelOK),velDiscret2,omegaDiscret);
        

        if size(cAngle1(indiceVelOK))~=size(omega1(indiceVelOK))
            cAngle1=cAngle1';
        end
        if size(cAngle2(indiceVelOK))~=size(omega2(indiceVelOK))
            cAngle2=cAngle2';
        end

		
		fct((indiceVelOK),1)  = ((omega2(indiceVelOK).*cAngle2(indiceVelOK)-omega1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),2)  = ((omega2(indiceVelOK).*cAngle2(indiceVelOK)-omegaLow1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),3)  = ((omega2(indiceVelOK).*cAngle2(indiceVelOK)-omegaHigh1(indiceVelOK).*cAngle1(indiceVelOK)).^2)';
		fct((indiceVelOK),4)  = ((omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)-omega1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),5)  = ((omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)-omegaLow1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),6)  = ((omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)-omegaHigh1(indiceVelOK).*cAngle1(indiceVelOK)).^2)';  
		fct((indiceVelOK),7)  = ((omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)-omega1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),8)  = ((omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)-omegaLow1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),9)  = ((omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)-omegaHigh1(indiceVelOK).*cAngle1(indiceVelOK)).^2)';
		if (length(indiceVelOK)<round(0.9*length(velData1)))  %  90% of the vel must be under the line
            fct(indiceVelNo,:)=c*(length(velData1)-length(indiceVelOK))^2;
		else
            fct(indiceVelNo,:)=1e-20;
		end
		
		
		fct=min(fct');
		
		fct=sum(fct);
		
    case 'fixed'
		c=0.1;
		
		
			c=0.1;
		
		
		% discretization
		domega=0.001;
		omegaDiscret=[0:domega:10];
        
        p1=(x(1)-((delta*x(2)*x(1))/vUnload))/x(2)+1;
        p2=(x(3)-((delta*x(4)*x(3))/vUnload))/x(4)+1;
        x=[x(1:2) p1 x(3:4) p2];
            
		
		velDiscret1=brShrinkDirect(x(1:3),omegaDiscret,delta);
		velDiscret2=brShrinkDirect(x(4:6),omegaDiscret,delta);
		velMax1=max(velDiscret1);
		velMin1=min(velDiscret1);
		velMax2=max(velDiscret2);
		velMin2=min(velDiscret2);
		
		
		%looking for the shrink velocties in the good range
		
		
		indiceVelOK=find(abs(velData1)>=velMin1 & abs(velData1)<=velMax1 & abs(velData2)>=velMin2 & abs(velData2)<=velMax2 & velData1<0 & velData2<0 ); 
		indiceVelNo=find(abs(velData1)<velMin1  | abs(velData1)>velMax1  | abs(velData2)<velMin2  | abs(velData2)>velMax2 |  velData1>0 | velData2>0 );
		
              
		
		
		velData1=abs(velData1);
		velData2=abs(velData2);
		
		% load calculation 
		
		[omega1(indiceVelOK),omegaHigh1(indiceVelOK),omegaLow1(indiceVelOK)] = brInvShrinkSpline(velData1(indiceVelOK),velDiscret1,omegaDiscret);
		[omega2(indiceVelOK),omegaHigh2(indiceVelOK),omegaLow2(indiceVelOK)] = brInvShrinkSpline(velData2(indiceVelOK),velDiscret2,omegaDiscret);
  
        
        if size(cAngle1(indiceVelOK))~=size(omega1(indiceVelOK))
            cAngle1=cAngle1';
        end
        if size(cAngle2(indiceVelOK))~=size(omega2(indiceVelOK))
            cAngle2=cAngle2';
        end

		
		fct((indiceVelOK),1)  = ((omega2(indiceVelOK).*cAngle2(indiceVelOK)-omega1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),2)  = ((omega2(indiceVelOK).*cAngle2(indiceVelOK)-omegaLow1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),3)  = ((omega2(indiceVelOK).*cAngle2(indiceVelOK)-omegaHigh1(indiceVelOK).*cAngle1(indiceVelOK)).^2)';
		fct((indiceVelOK),4)  = ((omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)-omega1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),5)  = ((omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)-omegaLow1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),6)  = ((omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)-omegaHigh1(indiceVelOK).*cAngle1(indiceVelOK)).^2)';  
		fct((indiceVelOK),7)  = ((omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)-omega1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),8)  = ((omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)-omegaLow1(indiceVelOK).*cAngle1(indiceVelOK)).^2)' ;
		fct((indiceVelOK),9)  = ((omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)-omegaHigh1(indiceVelOK).*cAngle1(indiceVelOK)).^2)';
		if (length(indiceVelOK)<round(0.9*length(velData1)))  %  90% of the vel must be under the line
            fct(indiceVelNo,:)=c*(length(velData1)-length(indiceVelOK))^2;
		else
            fct(indiceVelNo,:)=1e-20;
		end
		
		
		fct=min(fct');
		
		fct=sum(fct);
        
end
