% CC tracking 
% Liya 
% Oct 22, 2012

color_array= [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1;0 1 1 ;rand(20,3)];
iFrame = 1;
currentImg = double(MD.channels_(1).loadImage(iFrame));
figure(1);imagescc(currentImg);
finished_marking_flag=0;
nCell=0;
while(finished_marking_flag==0)
    h = imrect;
    position = wait(h);
    if(position(3)<1)
        % if the user input just 1 points, means no more cutting is needed
        finished_marking_flag=1;
    else
        nCell = nCell+1;
        position(3:4) = 2*round(position(3:4)/2)+1;
        position_array{iFrame}(nCell,1:4) = round(position);
    end
end

dmax = [20 20];
subpix = 'none';
d0=[0 0];
mkdir([MD.outputDirectory_,'/CC_tracking/']); 
img_width = size(currentImg,2);
img_height = size(currentImg,1);
pad_xy = [img_width/2 img_height/2];

for iFrame = 1 : MD.nFrames_;   
    
    previoucurrentImg = currentImg;
    currentImg = double(MD.channels_(1).loadImage(iFrame));
    
    h3 = figure(3);hold off; imagescc(currentImg);hold on;
    
    for iCell = 1 : nCell
        if iFrame >1
            previous_position = position_array{iFrame-1}(iCell,1:4);
            
            tPos = previous_position(1:2)+(previous_position(3:4)+1)/2;
            tDim = previous_position(3:4);
            
            displacement = ccbased_track(pad_boundary(previoucurrentImg),...
                [tPos(1) tPos(2)]+pad_xy,[tDim(1) tDim(2)],pad_boundary(currentImg),dmax,subpix,d0);
            
            current_position(3:4) = previous_position(3:4);
            current_position(1:2) = previous_position(1:2)+displacement;
            
            if current_position(1)<=0
                current_position(3) = 2*(round((current_position(3) + current_position(1) -2)/2)) + 1;
                current_position(1) = 1;
            end
            
            if current_position(2)<=0
                current_position(4) = 2*(round((current_position(4) + current_position(2) -2)/2)) + 1;
                current_position(2) = 1;
            end
            
            if current_position(1)+ current_position(3) > size(currentImg,2)-1
                current_position(3) = 2*(round((size(currentImg,2)-2 - current_position(1))/2)) + 1;
            end
            
            if current_position(2)+ current_position(4) > size(currentImg,1)-1
                current_position(4) = 2*(round((size(currentImg,1)-2 - current_position(2))/2)) + 1;
            end
            
            tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
         
            position_array{iFrame}(iCell,1:4) = current_position;
        else
            current_position = position_array{iFrame}(iCell,1:4);
            tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
        end
        plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        X = [current_position(1);...
            current_position(1)+ current_position(3);...
            current_position(1)+current_position(3);...
            current_position(1);...
            current_position(1);];
        
        Y = [current_position(2);...
            current_position(2);...
            current_position(2)+ current_position(4);...
            current_position(2)+current_position(4);...
            current_position(2);];
        plot(X,Y,'color',color_array(iCell,:));
        
    end
    title(['Tracking Frame ',num2str(iFrame)]);
    saveas(h3,[MD.outputDirectory_,'/CC_tracking/tracking_',num2str(iFrame),'.tif']);
end

save('position_array','position_array');

currentImg = double(MD.channels_(1).loadImage(1));

h3 = figure(3);hold off; imagescc(currentImg);hold on;
   
for iFrame = 1 : MD.nFrames_    
    
    for iCell = 1 : nCell
        
        current_position = position_array{iFrame}(iCell,1:4);
        tracked_pos = current_position(1:2)+(current_position(3:4)+1)/2;
        
        plot(tracked_pos(1),tracked_pos(2),'.','color',color_array(iCell,:));
        X = [current_position(1);...
            current_position(1)+ current_position(3);...
            current_position(1)+current_position(3);...
            current_position(1);...
            current_position(1);];
        
        Y = [current_position(2);...
            current_position(2);...
            current_position(2)+ current_position(4);...
            current_position(2)+current_position(4);...
            current_position(2);];
        plot(X,Y,'color',color_array(iCell,:));
        if iFrame >1 
            previous_position = position_array{iFrame-1}(iCell,1:4);
            tracked_old_pos = previous_position(1:2)+(previous_position(3:4)+1)/2;
            quiver(tracked_old_pos(1), tracked_old_pos(2), tracked_pos(1)-tracked_old_pos(1), tracked_pos(2)-tracked_old_pos(2),'color',color_array(iCell,:));
        
        end
    end
    

end

saveas(h3,[MD.outputDirectory_,'/CC_tracking/whole_tracking.tif']);