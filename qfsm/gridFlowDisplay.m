function gridFlowDisplay(dataPath,roidim,framerate,res,imgPathName,nFrames)
% gridFlowDisplay(dataPath,roidim)
% This function displays flow vectors at grids on top of image.
% example: gridFlowDisplay('/home/sh268/Documents/Alexis/MyosinII_Flows/Second_Analysis/Jasp/Jasp2',
%                         [w,h],5,72,'./data/Jasp2sdc_',30)
% input:
%       dataPath:   path to the folder
%       roidim =    [w,d] width (w) and height (h) of the roi dimension for
%                   cropping and quantifying the flow
%       framerate:  sec/frame
%       res:        image resolution, nm/pixel
%       imgPathName:path and name of myosin image 
%       nFrames:    total frame number
% output:
%       images of flow vectors (.fig):  stored in './roinum/flowimg/'
%       flow magnitudes in the roi:     maxflow (after outlier) mag per
%                                       each frame, stored in
%                                       './roinum/maxflow'
% Sangyoon Han Feb 2013
cd(dataPath)
flow = []; % flow initialization
load([dataPath '/corr/flow/flow000.mat'])
nNodes = length(flow);
% nFrames = 28;
% res = 72; % nm/pix
dt = framerate; % sec/frame
width = roidim(1);
height = roidim(2);
avgFlowSpeed = zeros(nFrames,1);
maxFlowSpeed = zeros(nFrames,1);
% imgPathName = './data/Jasp2sdc_';
% imgPathName = './data/stackMRLC_';
bshow = true;%false;
clear flow
%% loop
if bshow
    h = figure;
end
for k=0:nFrames-1
    % setup
    kformat = ['%.' '3' 'd'];
    flow = load(strcat(dataPath,'/corr/flow/flow',num2str(k,kformat),'.mat'),'flow');
    flow = flow.flow;
    stackImg = imread(strcat(imgPathName,num2str(k,kformat),'.tif'));

    if bshow
        % vector filtering
        % Masking with the myosin image
        grayThres = graythresh(stackImg);
        bwstackImg = im2bw(stackImg,grayThres);
        bwstackImg = bwmorph(bwstackImg,'clean',Inf);
        outerMask = bwmorph(bwstackImg,'close',1);
        bwstackImg = bwmorph(bwstackImg,'dilate',15);
        bwstackImg = bwmorph(bwstackImg,'erode',23);
        innerMask = imcomplement(bwstackImg);
        bwstackImg = innerMask & outerMask;
        % filtering out the flow vectors outside the mask
        idx = true(1,numel(flow(:,1)))'; %index for filtering

        outsideIdx = arrayfun(@(i,j) maskVectors(i,j,bwstackImg),flow(:,2),flow(:,1));
        % outlier detection
        dispMat = [flow(:,2),flow(:,1),flow(:,4)-flow(:,2),flow(:,3)-flow(:,1)];
        outlierThreshold = 2;
        outlierIndex = detectVectorFieldOutliers(dispMat,outlierThreshold,1);
        idx(outlierIndex) = false;
        
        idx = idx & outsideIdx;
        flow_filtered = flow(idx,:);
        flow(~idx,3:4) = NaN(sum(~idx),2);
                
        % showing    
        hImg = imshow(stackImg,[]);
        hold on
        mag = 100/dt;
        hflow = quiver(flow_filtered(:,2),flow_filtered(:,1),mag*(flow_filtered(:,4)-flow_filtered(:,2)),mag*(flow_filtered(:,3)-flow_filtered(:,1)),0,'r','Linewidth',0.5);
        if k==0
            % draw line for rotation
%             uiwait(msgbox('Draw a line for roi to set rotation angle. Click OK to continue.','mage Rotation','modal'));
%             [rot_y1 rot_x1]=ginput(1);
%             plot(rot_y1,rot_x1,'b+');
%             [rot_y2 rot_x2]=ginput(1);
%             plot(rot_y2,rot_x2,'b+');
% 
%             rot_angle = -atan2((rot_x1-rot_x2),(rot_y2-rot_y1))*180/pi;
%             Xr=imrotate(X,rot_angle,'bilinear','crop');

            % cropping the roi
            display('Please draw a rectangle of roi')
            hrect = imrect;
            posRect = wait(hrect); %[xmin ymin w h]
            delete(hrect)
            centerx = posRect(1) + posRect(3)/2;
            centery = posRect(2) + posRect(4)/2;
            xmin = centerx - width/2;
            ymin = centery - height/2;
            xmax = centerx + width/2;
            ymax = centery + height/2;
            roinum = input('ROI number for this cell?');
        end

        % zoom in
        Coord = get(h,'Position');
        if width>=height
            set(h,'Position',[Coord(1),Coord(2),720,276])
        else
            set(h,'Position',[Coord(1),Coord(2),276,720])
        end            
        xlim([xmin xmin+width]);
        ylim([ymin ymin+height]);

        delete(hflow)
        hflow = quiver(flow_filtered(:,2),flow_filtered(:,1),mag*(flow_filtered(:,4)-flow_filtered(:,2)),mag*(flow_filtered(:,3)-flow_filtered(:,1)),0,'g','Linewidth',1);

        % save flow magnitudes:
        dispMat = [flow(:,2),flow(:,1),flow(:,4)-flow(:,2),flow(:,3)-flow(:,1)];
        flow_mag = (dispMat(:,3).^2+dispMat(:,4).^2).^0.5*res/dt;
        flowMagAll(:,k+1) = flow_mag;

        % second roi
        if k==0
            % cropping the roi
            display('Please draw an arc for cortical actin boundary clockwise. Finish with right click after a final point')
%             delete(hrect);
            hpoly = impoly;
            posPoly = wait(hpoly); 
            delete(hpoly);
        end
        
        % I need to filter out vectors that point out from center of a cell
        % for which I need to compare orientation-based filter.
        % finding line vectors perpendicular to the drawn lines
        vecArcs = zeros(2,length(posPoly)-1);
        for ii=1:length(posPoly)-1
            dline = posPoly(ii:ii+1,:);
            vecArcs(:,ii)=[0 -1;1 0]*[dline(2,1)-dline(1,1); dline(2,2)-dline(1,2)]; %[-pi pi]
            vecArcs(:,ii) = vecArcs(:,ii)/norm(vecArcs(:,ii));% normalize
        end
        % determining the range of angles
        % make angles continuous
        vec_middle = [mean(vecArcs(1,:));mean(vecArcs(2,:))];
        % in when angles are in between ang_middle-pi/4 and ang_middle+pi/4
        in  = inpolygon(flow_filtered(:,2),flow_filtered(:,1),[xmin xmax xmax xmin],[ymin ymin ymax ymax]);
        flowROI = flow_filtered(in,:);
        vecROI = [flowROI(:,4)-flowROI(:,2),flowROI(:,3)-flowROI(:,1)];

        inAngIdx = arrayfun(@(i,j) angleFilter(i,j,vec_middle),vecROI(:,1),vecROI(:,2));

        % comparing angles of vectors from this range
        flowFinal = flowROI(inAngIdx,:);
        % create mask that is outside of posPoly
%         in  = inpolygon(flow_filtered(:,2),flow_filtered(:,1),posPoly(:,1),posPoly(:,2));
%         flowROI = flow_filtered(in,:);
        delete(hflow)
        hflow = quiver(flowFinal(:,2),flowFinal(:,1),mag*(flowFinal(:,4)-flowFinal(:,2)),mag*(flowFinal(:,3)-flowFinal(:,1)),0,'g','Linewidth',1);
        % scale vector
        if width>=height
            hlscale = line([xmax-10 xmax-10+500/res],[ymax-2 ymax-2],'Color',[1 1 1],'Linewidth',2);
            hfscale = quiver(xmin+2,ymax-2,mag*500/res/dt,0,0,'w','Linewidth',0.5); % 500 nm/sec
%             ht = text(xmin+1+mag*500/res/dt/3,ymax-1,'10 nm/sec');
        else
            hlscale = line([xmax-2 xmax-2],[ymax-10 ymax-10+500/res],'Color',[1 1 1],'Linewidth',2);
            hfscale = quiver(xmax-2,ymin+2,0,mag*500/res/dt,0,'w','Linewidth',1); % 500 nm/sec
%             ht = text(xmin+1+mag*500/res/dt/3,ymax-1,'10 nm/sec');
        end
%         set(ht,'Color','g')
%         set(ht,'FontName','Arial')
%         set(ht,'FontSize',7)
        
        ROIpath = ['./results/roi' num2str(roinum)];
        imgflowPath = [ROIpath '/imgflow/'];
        flowdataPath = [ROIpath '/flowdata/'];
        
        if k==0 && (~exist(imgflowPath,'dir') || ~exist(flowdataPath,'dir'))
            mkdir(imgflowPath);
            mkdir(flowdataPath);
        end
%         set(h,'Position',[Coord(1),Coord(2),720,276])

        hgexport(h,strcat(imgflowPath,'flowTif',num2str(k,kformat)),hgexport('factorystyle'),'Format','tiff')
        hgsave(h,strcat(imgflowPath,'flowFig',num2str(k,kformat)),'-v7.3')
        save(strcat(flowdataPath,'flow_ROI',num2str(k,kformat),'.mat'),'flowROI')
        delete(hflow)
        delete(hImg)
        delete(hfscale) 
        delete(hlscale)
%         delete(ht)
    end
    % collecting flow magnitudes
    % saving each magnitude 
%       flow magnitudes in the roi:     maxflow (after outlier) mag per
%                                       each frame, stored in
%                                       './roinum/maxflow'
    flow_mag = ((flowFinal(:,4)-flowFinal(:,2)).^2+(flowFinal(:,3)-flowFinal(:,1)).^2).^.5;
    avgSpeed = nanmean(flow_mag)*res/dt;
    maxSpeed = nanmax(flow_mag)*res/dt;
    avgFlowSpeed(k+1) = avgSpeed;
    maxFlowSpeed(k+1) = maxSpeed;
    save(strcat(flowdataPath,'flow_mag',num2str(k,kformat),'.mat'),'flow_mag');
    save(strcat(flowdataPath,'flow',num2str(k,kformat),'.mat'),'flowFinal');
%     save(strcat(flowdataPath,'avgOuterFlowSpeed',num2str(k,kformat),'.mat'),'avgSpeed');
%     % 3D matrix making
%     fmagmat(:,:,k+1) = [flow(:,2),flow(:,1),flow_mag];
end
save(strcat(flowdataPath,'flowMagAll.mat'),'flowMagAll');

% %% average speed for outer flow
% nanmean(flow_mag)*res/dt
% nanstd(flow_mag)*res/dt
%% max flow mag for entire flow
% maxFlow = squeeze(nanmax(fmagmat(:,3,:),[],1));
display(['mean of avg flow is ' num2str(nanmean(avgFlowSpeed)) ' nm/sec.'])
display(['standard deviation is ' num2str(nanstd(avgFlowSpeed)) ' nm/sec.'])
display(['mean of max flow is ' num2str(nanmean(maxFlowSpeed)) ' nm/sec.'])
display(['standard deviation is ' num2str(nanstd(maxFlowSpeed)) ' nm/sec.'])
save(strcat(flowdataPath,'avgOuterFlowSpeed.mat'),'avgFlowSpeed');
save(strcat(flowdataPath,'maxOuterFlowSpeed.mat'),'maxFlowSpeed');
return
%% min flow mag for entire flow
minFlow = squeeze(nanmin(fmagmat(:,3,:),[],1));
display(['mean of min flow is ' num2str(nanmean(minFlow)*res/dt) ' nm/sec.'])
display(['standard deviation is ' num2str(nanstd(minFlow)*res/dt) ' nm/sec.'])
%% avg flow mag for entire flow
avgFlow = squeeze(nanmean(fmagmat(:,3,:),1));
display(['mean of min flow is ' num2str(nanmean(avgFlow)*res/dt) ' nm/sec.'])
display(['standard deviation is ' num2str(nanstd(avgFlow)*res/dt) ' nm/sec.'])
%% Heatmap of average flow speed
% load all flow vector in to 3D array
imshow(stackImg,[])
hold on
% averaging the flow mag
avgMag = nanmean(fmagmat,3);
[X,Y]=meshgrid(unique(avgMag(:,1)),unique(avgMag(:,2)));
Z = zeros(size(X));
[nrows,ncols]=size(X);
for jj=1:ncols
    for ii=1:nrows
        posx = X(ii,jj);
        posy = Y(ii,jj);
        % see if x,y position is in avgMag(:,1) and avgMag(:,2)
        indposx = avgMag(:,1)==posx;
        indposy = avgMag(:,2)==posy;
        ind = indposx & indposy;
        if sum(ind)==1
            Z(ii,jj) = avgMag(find(ind),3);
        elseif sum(ind)>1
            error('there are more than one flow at the same position');
        else
            Z(ii,jj) = NaN;
        end
    end
end
surf(X,Y,Z);view(0,90)
shading interp
colormap('jet')
axis equal;
view([0 0 1])
caxis auto
zlim([0 5])

% mesh(avgMag(:,1),avgMag(:,2),avgMag(:,3))
title('Average flow heatmap') 
set(gca, 'DataAspectRatio', [1,1,1]);
%% zooming
xlim([633 728]);
ylim([283 356]);
%% Plotting the magnitudes of the flow
% average flow magnitude
t = 1:nFrames;
nValidNodes = 0;
avgFlow = zeros(1,nFrames);
for ii=1:nNodes
    %only numerical arrays are stored
    if ~any(~flow_index(ii,:))
        nValidNodes = nValidNodes+1;
        flow_mag_val(nValidNodes,:) = flow_mag(ii,:);
        avgFlow = avgFlow + flow_mag_val(nValidNodes,:);
        plot(t,flow_mag_val(nValidNodes,:))
        hold all
    end
end
avgFlow = avgFlow/nValidNodes;
plot(t,avgFlow,'k','Linewidth',4)
