function [fct]=brBothDiffOptim(x,velData1,velData2,cAngle1,cAngle2,delta,vUnload,pFixe);

% Looking which mt is growing

if length(find(velData1<0))==length(velData1) & length(find(velData2>0))==length(velData2)
    state='SG';
elseif length(find(velData1>0))==length(velData1) & length(find(velData2<0))==length(velData2)
    state='GS';
else
    error('You do not input velocity of growing/shrinkage or shrinkage/growing case');
end

c=0.1;

switch pFixe
    
    case 'free'
		switch state
            
            case 'SG'
                
                
                
                % discretization
				domega=0.001;
				omegaDiscret=[0:domega:10];
				
				velDiscret1=brShrinkDirect(x(1:3),omegaDiscret,delta);
				velMax1=max(velDiscret1);
				velMin1=min(velDiscret1);
		
                
                % range control
                
                indiceVelOK=find(abs(velData1)>=velMin1 & abs(velData1)<=velMax1 & velData1<0 & velData2>0 ); 
                indiceVelNo=find(abs(velData1)<velMin1  | abs(velData1)>velMax1  | velData1>0 | velData2<0 );
                
                % Load calculation
                velData1=abs(velData1);
                [omega1(indiceVelOK),omegaHigh1(indiceVelOK),omegaLow1(indiceVelOK)] = brInvShrinkSpline(velData1(indiceVelOK),velDiscret1,omegaDiscret);
                omega2(indiceVelOK)                                                  = peskinInv(velData2(indiceVelOK),x(4:5),delta);
                
                if size(cAngle1(indiceVelOK))~=size(omega1(indiceVelOK))
                    cAngle1=cAngle1';
                end
                if size(cAngle2(indiceVelOK))~=size(omega2(indiceVelOK))
                    cAngle2=cAngle2';
                end

                fct((indiceVelOK),1)  = ((cAngle2(indiceVelOK).*omega2(indiceVelOK)-cAngle1(indiceVelOK).*omega1(indiceVelOK)).^2)' ;
				fct((indiceVelOK),2)  = ((cAngle2(indiceVelOK).*omega2(indiceVelOK)-cAngle1(indiceVelOK).*omegaLow1(indiceVelOK)).^2)' ;
				fct((indiceVelOK),3)  = ((cAngle2(indiceVelOK).*omega2(indiceVelOK)-cAngle1(indiceVelOK).*omegaHigh1(indiceVelOK)).^2)';
                if (length(indiceVelOK)<round(0.9*length(velData1)))  %  90% of the vel must be under the line
                    fct(indiceVelNo,:)=c*(length(velData1)-length(indiceVelOK))^2;
				else
                    fct(indiceVelNo,:)=1e-20;
				end
                
                fct=min(fct');
                fct=sum(fct);
      
        case 'GS'
                    % discretization
			domega=0.001;
			omegaDiscret=[0:domega:10];
			
			velDiscret2=brShrinkDirect(x(3:5),omegaDiscret,delta);
			velMax2=max(velDiscret2);
			velMin2=min(velDiscret2);
	
            
            % range control
            
            indiceVelOK=find(abs(velData2)>=velMin2 & abs(velData2)<=velMax2 & velData2<0 & velData1>0 ); 
            indiceVelNo=find(abs(velData2)<velMin2  | abs(velData2)>velMax2  | velData2>0 | velData1<0 );
            
            velData2=abs(velData2);
            % Load calculation
            
            [omega2(indiceVelOK),omegaHigh2(indiceVelOK),omegaLow2(indiceVelOK)] = brInvShrinkSpline(velData2(indiceVelOK),velDiscret2,omegaDiscret);
            omega1(indiceVelOK)                                                  = peskinInv(velData1(indiceVelOK),x(1:2),delta);
            if size(cAngle1(indiceVelOK))~=size(omega1(indiceVelOK))
                cAngle1=cAngle1';
            end
            if size(cAngle2(indiceVelOK))~=size(omega2(indiceVelOK))
                cAngle2=cAngle2';
            end                                                  
            fct((indiceVelOK),1)  = ((omega1(indiceVelOK).*cAngle1(indiceVelOK)-omega2(indiceVelOK).*cAngle2(indiceVelOK)).^2)' ;
			fct((indiceVelOK),2)  = ((omega1(indiceVelOK).*cAngle1(indiceVelOK)-omegaLow2(indiceVelOK).*cAngle2(indiceVelOK)).^2)' ;
			fct((indiceVelOK),3)  = ((omega1(indiceVelOK).*cAngle1(indiceVelOK)-omegaHigh2(indiceVelOK).*cAngle2(indiceVelOK)).^2)';
            if (length(indiceVelOK)<round(0.9*length(velData1)))  %  90% of the vel must be under the line
                fct(indiceVelNo,:)=c*(length(velData1)-length(indiceVelOK))^2;
			else
                fct(indiceVelNo,:)=1e-20;
			end
            
            fct=min(fct');
            fct=sum(fct);
        end
                
    case 'fixed'
		switch state
            
            case 'SG'
                
                % discretization
				domega=0.001;
				omegaDiscret=[0:domega:10];
                
                p1=(x(1)-((delta*x(2)*x(1))/vUnload))/x(2)+1;
				
                x=[x(1:2) p1 x(3:4)];
                
				velDiscret1=brShrinkDirect(x(1:3),omegaDiscret,delta);
				velMax1=max(velDiscret1);
				velMin1=min(velDiscret1);
		
                
                % range control
                
                indiceVelOK=find(abs(velData1)>=velMin1 & abs(velData1)<=velMax1 & velData1<0 & velData2>0 ); 
                indiceVelNo=find(abs(velData1)<velMin1  | abs(velData1)>velMax1  | velData1>0 | velData2<0 );
                
                % Load calculation
                velData1=abs(velData1);
                [omega1(indiceVelOK),omegaHigh1(indiceVelOK),omegaLow1(indiceVelOK)] = brInvShrinkSpline(velData1(indiceVelOK),velDiscret1,omegaDiscret);
                omega2(indiceVelOK)                                                  = peskinInv(velData2(indiceVelOK),x(4:5),delta);
                
                if size(cAngle1(indiceVelOK))~=size(omega1(indiceVelOK))
                    cAngle1=cAngle1';
                end
                if size(cAngle2(indiceVelOK))~=size(omega2(indiceVelOK))
                    cAngle2=cAngle2';
                end

                fct((indiceVelOK),1)  = ((cAngle2(indiceVelOK).*omega2(indiceVelOK)-cAngle1(indiceVelOK).*omega1(indiceVelOK)).^2)' ;
				fct((indiceVelOK),2)  = ((cAngle2(indiceVelOK).*omega2(indiceVelOK)-cAngle1(indiceVelOK).*omegaLow1(indiceVelOK)).^2)' ;
				fct((indiceVelOK),3)  = ((cAngle2(indiceVelOK).*omega2(indiceVelOK)-cAngle1(indiceVelOK).*omegaHigh1(indiceVelOK)).^2)';
                if (length(indiceVelOK)<round(0.9*length(velData1)))  %  90% of the vel must be under the line
                    fct(indiceVelNo,:)=c*(length(velData1)-length(indiceVelOK))^2;
				else
                    fct(indiceVelNo,:)=1e-20;
				end
                
                fct=min(fct');
                fct=sum(fct);
                
            case 'GS'
                        % discretization
				domega=0.001;
				omegaDiscret=[0:domega:10];
                p2=(x(3)-((delta*x(4)*x(3))/vUnload))/x(4)+1;
                
                x=[x(1:2)  x(3:4) p2];
				
				velDiscret2=brShrinkDirect(x(3:5),omegaDiscret,delta);
				velMax2=max(velDiscret2);
				velMin2=min(velDiscret2);
		
                
                % range control
                
                indiceVelOK=find(abs(velData2)>=velMin2 & abs(velData2)<=velMax2 & velData2<0 & velData1>0 ); 
                indiceVelNo=find(abs(velData2)<velMin2  | abs(velData2)>velMax2  | velData2>0 | velData1<0 );
                
                velData2=abs(velData2);
                % Load calculation
                
                [omega2(indiceVelOK),omegaHigh2(indiceVelOK),omegaLow2(indiceVelOK)] = brInvShrinkSpline(velData2(indiceVelOK),velDiscret2,omegaDiscret);
                omega1(indiceVelOK)                                                  = peskinInv(velData1(indiceVelOK),x(1:2),delta);
                
                if size(cAngle1(indiceVelOK))~=size(omega1(indiceVelOK))
                    cAngle1=cAngle1';
                end
                if size(cAngle2(indiceVelOK))~=size(omega2(indiceVelOK))
                    cAngle2=cAngle2';
                end

                                                      
                fct((indiceVelOK),1)  = ((cAngle1(indiceVelOK).*omega1(indiceVelOK)-cAngle2(indiceVelOK).*omega2(indiceVelOK)).^2)' ;
				fct((indiceVelOK),2)  = ((cAngle1(indiceVelOK).*omega1(indiceVelOK)-cAngle2(indiceVelOK).*omegaLow2(indiceVelOK)).^2)' ;
				fct((indiceVelOK),3)  = ((cAngle1(indiceVelOK).*omega1(indiceVelOK)-cAngle2(indiceVelOK).*omegaHigh2(indiceVelOK)).^2)';
                if (length(indiceVelOK)<round(0.9*length(velData1)))  %  90% of the vel must be under the line
                    fct(indiceVelNo,:)=c*(length(velData1)-length(indiceVelOK))^2;
				else
                    fct(indiceVelNo,:)=1e-20;
				end
                
                fct=min(fct');
                fct=sum(fct);
		end
end
