function [coor]= takenkick(info,coor)
% takenkick displays an image with coordinates ploted into it. It allows to
%           manually add and delete coordinates
%
% SYNOPSIS       [coor]= takenkick(info,coor)
%
% INPUT          info : image
%                coor : list of coordinates
%
%
% OUTPUT         coor : edited list of coordinates
%                             
%
% DEPENDENCIES   takenkick uses {nothing}
%                                  
%                takenkick is used by { ptInitializeJob }
%
% Colin Glass, Feb 04

% The goal of this programm is to manually add or substract coordinates.
% So a picture will pop up whilst running the programm and the user clicks on
% marks missed or faulty coordinates


taken = figure, imshow(info,[]), title('Left-klick on the unmarked nuclei, right-klick on the falsely marked nuclei /// Press ENTER twice when finished.');
hold on;

plot(coor(:,1),coor(:,2),'r.');

done=0;
found=[];
addon=[0,0];
counter=0;

while done==0
        %the user can click on cells. Everytime he presses ENTER, his
        %changes will become visible.
        %He can continue to mark cells as long as long as he marks at 
        %least one cell between pressing ENTER twice.
        %Then he will be asked if this is really correct. Leftclick
        %indicates yes, stopping the programm. Rightclick will let him
        %continue clicking aroung like a madman.
       
		counter=counter+1;
		hold on;
		xx=[];
		yy=[];
		leftright=[];
		[xx,yy,leftright] = ginput;
		
		%if no Input - finish
		if isempty(xx)
            
                 %wipe out the cells marked as faulty
                 if isempty (found)==0
                       coor(found,:)=[];
                 end
                 
                 %add the desired extra cells
                 addon(1,:)=[];
                 coor=cat(1,coor,addon);
                 
                 %ask the user if this is really correct
                 onoff=[];
                 onoff=0;
                 figure, imshow(info,[]), title('Is this correct? if yes: leftklick  -- if no: rightklick    then press ENTER')
                 hold on;
                 plot(coor(:,1),coor(:,2),'r.');
                 [nothing,nothing2,onoff] = ginput;
                 hold off
                 
                 close
                 
                 %onoff is one if the user did a left click. In this case
                 %the programm will stop.
                 if isempty (onoff)
                     done=1;
                     break
                     
                 elseif onoff(1)==1
                       done=1;
                       break
                 end
                 
                 %if he didn't, the programm goes on offering him to change
                 %as he wishes. Since the addon's and the delates up to
                 %this moment are already integrated we have to
                 %reinitialize these two variables
                 addon=[0,0];
                 found=[];
                 
                 figure, imshow(info,[]), title('LEFT klick on the UNMARKED(!!!!!) nucloi /// RIGHT klick on (or close to) the FALSELY marked nucs ///  press ENTER when finished /// repeats as long as you klick on something, between to sucsessive presses on ENTER...');
                 hold on;
                 plot(coor(:,1),coor(:,2),'r.');
                 [xx,yy,leftright] = ginput;
                 close
                 
		end
        %as long as there's Input - continue
        
		left=[];
		right=[];
		
        %distinguish left from rightclicks
		left=find(leftright==1);
		right=find(leftright==3);
         
		if isempty(left)==0
                 %save the left clicks (cells to be added) and plot them
                 %into the picture
                 cooradd=[0,0];
                 xin=round(xx(left));
                 yin=round(yy(left));
                 plot(xin,yin,'r.');
                 cooradd=cat(2,xin,yin);
                 addon=cat(1,addon,cooradd);
		end
		
		
		if isempty(right)==0
                  %find the cells that should be removed (nearest to a rightclick) and plot a yellow spot
                  %on top of them for visualization
                  xout=round(xx(right));
                  yout=round(yy(right));
			      
                  %finding...
                  for i=1:length(xout)
                        [noth,found(end+1,1)]=min(abs(coor(:,1)-xout(i)) + abs(coor(:,2)-yout(i)));
                  end
                  
                  coorout=[0,0];
                  coorout=coor(found,:);
                  
                  hold on;
                  plot(coorout(:,1),coorout(:,2),'y.');
                  hold off;
		end

		
		hold off;

end
close(taken);
