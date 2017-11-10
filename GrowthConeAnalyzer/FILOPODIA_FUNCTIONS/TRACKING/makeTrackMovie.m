function [outputTracking] = makeTrackMovie(tracksCoordAmpLink,trackedFeatureIndex,analInfo, imDir,saveDir,timeInt,pixSize_um)

% should have three different options : make movie, get timeseries, vel or
% length 

if isempty(imDir)
    imDir = uigetdir(pwd,'select a directory where the images are stored');
end 

if isempty(saveDir)
   saveDir  = uigetdir(pwd,'select a directory to save the results');
end 
    

% collect images and initiate
[listOfImages] = searchFiles('.tif',[],imDir,0);
nImTot = size(listOfImages,1);

%% sort images if not padded 
 for iFrame =  1:nImTot;
        imageName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
        
        %Call "getFilenameBody" from common dir to split filename into path,
        %body, no, and ext. Put path, body and ext into appropriate cell/number vector for that
        %frame
        [path body no extsort ] = getFilenameBody(imageName);
        
        
        pathCell(iFrame) = cellstr(path);
        bodyCell(iFrame) = cellstr(body);
        extCell(iFrame) = cellstr(extsort);
        
        % output "no" is a string so convert to number to sort
        num(iFrame)  = str2double(no);
        
  end
%    Sort number vector numerically
    sortednum = sort(num);
    
    %Convert number vector to cell
    sortednum_cell = num2cell(sortednum);
    
    %Create Sorted Image List
    sortedImages = [pathCell', bodyCell', sortednum_cell', extCell'];

%%





%Configure Figure 
fileName = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
%fileName = [char(sortedImages(1,1)) filesep char(sortedImages(1,2)),...;
                  %num2str(sortednum(1)) char(sortedImages(1,4))];
img = double(imread(fileName)); 
img = -img;
img = [img img];
[ny nx] = size(img) ; 
zoom = 1; 
ext = '.png';

%% for now make sure to add extra filter field.. should technically do this 
   % right after fit
if ~isfield(analInfo(1).filoInfo,'Ext_filoPixCut'); 
   % prune the filo based on endpoint coord *not elegant but quick fix for
   % now!! 
   for iFrame = 1:numel(analInfo)
       % get all filo Info for the frame
          filoInfoCurrent = analInfo(iFrame).filoInfo; 
       
       for ifilo = 1: numel(filoInfoCurrent)
    
     
       endpoint = filoInfoCurrent(ifilo).Ext_endpointCoordFitPix; 
       filoPix = filoInfoCurrent(ifilo).Ext_pixIndices; 
       cutIdx  = find(filoPix==endpoint); 
       filoPixCut = filoPix(1:cutIdx); 
       % save the new field
       analInfo(iFrame).filoInfo(ifilo).Ext_filoPixCut = filoPixCut; 
       end 
   end 
end    

% filter tracks coords 
trackInfo= getTrackSEL(tracksCoordAmpLink);  

% get all x and yCoords (row = track, col = frame)
xCoordsPreFilt  = tracksCoordAmpLink(:,1:8:end); 
yCoordsPreFilt = tracksCoordAmpLink(:,2:8:end);
numTracks = size(tracksCoordAmpLink);


% filter all tracks less than 3 frames
tracksCoordAmpLink(trackInfo(:,3)<3,:) = []; 
  lengths = tracksCoordAmpLink(:,4:8:end);
  lengths = lengths.*pixSize_um;% in um
  vels = diff(lengths,1,2);
  vels=vels./timeInt*60; % in convert to um/min

%% MAKE COLORMAPS OF TIME SERIES : AVG VEL EACH TRACK
  % make average velocity map 
  velsMean = nanmean(vels,2); 
  
  avgVel = repmat(velsMean,1,length(vels(1,:))); 
  
  avgVel  = avgVel.*~isnan(vels); 
  maxSpeedAvg = prctile(velsMean(:),95); 
  minSpeedAvg = prctile(velsMean(:),5); 
  
  avgVel(avgVel==0)= NaN;  
  figure; 
  imagesc(avgVel)
  colorbar; 
  set(gca,'CLim',[minSpeedAvg,maxSpeedAvg]); 
  
  % Create xlabel
    xlabel({'Frame Number'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');

% Create ylabel
    ylabel({'Filopodia Track ID'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');
  
 saveas(gcf,[saveDir filesep 'meanVelocityTrackAll.eps'],'psc2'); 
 saveas(gcf,[saveDir filesep 'meanVelocityTrackAll.tif']);

close gcf 
%%  MAKE COLORMAPS OF TIME SERIES : FRAME TO FRAME VEL EACH TRACK
  maxSpeed = prctile(vels(:),95);
  minSpeed = prctile(vels(:),5);

  
    imagesc(vels);
     colorbar;
    set(gca,'CLim',[minSpeed ,maxSpeed]);
   
 
  
  % Create xlabel
    xlabel({'Frame Number'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');

% Create ylabel
    ylabel({'Filopodia Track ID'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');
  saveas(gcf,[saveDir filesep 'velocityAll.eps'],'psc2'); 
  saveas(gcf,[saveDir filesep 'velocityAll.tif']);
  close gcf
%%  MAKE COLORMAPS OF TIME SERIES : LENGTHS OF EACH TRACK

maxLength = prctile(lengths(:),95); 
minLength = prctile(lengths(:),5);

imagesc(lengths) 
colorbar; 
set(gca,'CLim',[minLength,maxLength]); 

  % Create xlabel
    xlabel({'Frame Number'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');

% Create ylabel
    ylabel({'Filopodia Track ID'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');
  saveas(gcf,[saveDir filesep 'lengthsAll.eps'],'psc2'); 
  saveas(gcf,[saveDir filesep 'lengthsAll.tif']);
  close gcf
  
  %%
% get lengths
numTracksFilt = size(tracksCoordAmpLink);
trackID = 1:numTracksFilt; 
xCoordsPostFilt = tracksCoordAmpLink(:,1:8:end); 
yCoordsPostFilt = tracksCoordAmpLink(:,2:8:end); 


%% filter further if desired 

getRetract = 0; 
if getRetract == 1
    trackInfo= getTrackSEL(tracksCoordAmpLink);
   tracksCoordAmpLink = tracksCoordAmpLink(velsMean<-0.5 &trackInfo(:,3) >10,:); 
    trackID = trackID(velsMean<-0.5&trackInfo(:,3)>10); 
    if ~isdir([saveDir filesep 'RetractMovie']) 
        mkdir([saveDir filesep 'RetractMovie']); 
    end 
    
    xCoordsPostFilt = tracksCoordAmpLink(:,1:8:end); 
    yCoordsPostFilt = tracksCoordAmpLink(:,2:8:end); 
    
    nTracksFilt = size(tracksCoordAmpLink); 
    
  
end 

getProtrusion = 0; 
if getProtrusion == 1
    trackInfo = getTrackSEL(tracksCoordAmpLink); 
    tracksCoordAmpLink = tracksCoordAmpLink(velsMean>0.5 &trackInfo(:,3) >7,:); 
    trackID = trackID(velsMean>0.5&trackInfo(:,3)>7); 
     xCoordsPostFilt = tracksCoordAmpLink(:,1:8:end); 
    yCoordsPostFilt = tracksCoordAmpLink(:,2:8:end); 
    
    nTracksFilt = size(tracksCoordAmpLink); 
end 

getStableFilos = 0 ;
if getStableFilos ==1
    
     trackInfo = getTrackSEL(tracksCoordAmpLink); 
    tracksCoordAmpLink = tracksCoordAmpLink(velsMean<0.5 & velsMean>-0.5&trackInfo(:,3) >7,:); 
    trackID = trackID(velsMean<0.5& velsMean>-0.5 &trackInfo(:,3)>7); 
     xCoordsPostFilt = tracksCoordAmpLink(:,1:8:end); 
    yCoordsPostFilt = tracksCoordAmpLink(:,2:8:end); 
    
    nTracksFilt = size(tracksCoordAmpLink); 
end 

for iFrame = 1:numel(analInfo)-1
   fileName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
  % fileName = [char(sortedImages(iFrame,1)) filesep char(sortedImages(iFrame,2)),...;
   %     num2str(sortednum(iFrame)) char(sortedImages(iFrame,4))];
    
    img = double(imread(fileName)); 
    img = -img;
    img = [img img];
    neuriteEdgeMask = analInfo(iFrame).masks.neuriteEdge; 
    roiYXNE = bwboundaries(neuriteEdgeMask); 
    
    
% reset handle
    h = figure('Visible','off','Position',[50 50 nx ny]) ;

iptsetpref('ImshowBorder','tight');

set(h, 'InvertHardcopy', 'off');

set(h, 'PaperUnits', 'Points');

set(h, 'PaperSize', [nx ny]);

set(h, 'PaperPosition', [0 0 nx ny]); % very important

%  set(h,'DefaultLineLineSmoothing','on');
%  set(h,'DefaultPatchLineSmoothing','on');

% Configure axes 
ha = axes('Position',[0 0 1 1]); 
fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];

   
    imshow(img,[]); 
hold on 
% plot all the detected filo ends for that frame
  %scatter(xCoordsPreFilt(:,iFrame),yCoordsPreFilt(:,iFrame),'filled','b');   
  % get all objects each frame
idxDetect  = trackedFeatureIndex(:,iFrame)';
idxDetect = idxDetect(idxDetect~=0); 
idxDetect = idxDetect(~isnan(idxDetect)); 
pixDetCurrFrame = vertcat(analInfo(iFrame).filoInfo(idxDetect).Ext_filoPixCut);
mask = zeros(size(img)); 
mask(pixDetCurrFrame) = 1; 
spy(mask,'k',5); 

% plot all the detected filo ends for that frame
%   scatter(xCoordsPreFilt(:,iFrame),yCoordsPreFilt(:,iFrame),'filled','b');
% trackedFeatureInfo
%  
%cMap=autumn(10);
 cMapLength=128; 
    cMap=mat2cell(jet(cMapLength),ones(cMapLength,1),3);
    mapper1=linspace(-8,8,cMapLength)';
    
   
    
for i = 1:length(trackID)
 % plot all the coords in each frame    
  % whatever just do it the stupid loopy way for now
    
   x = tracksCoordAmpLink(i,1:8:end); 
   % plot only tracks that are present in that frame 
   if ~isnan(x(iFrame))
     % calc mean velocity of that track   
     output.vel(i) = nanmean(vels(trackID(i),:));
     % map back to appropriate color
     
    D=createDistanceMatrix(output.vel(i),mapper1);
    [sD,idx]=sort(abs(D),2);
    idx=idx(:,1);
   
    colorTrack = cell2mat(cMap(idx,:)); 
     
  
      
   % plot the track path 
   plot(tracksCoordAmpLink(i,1:8:end),tracksCoordAmpLink(i,2:8:end),'LineStyle','-','Linewidth',1,'MarkerSize',5,'Color',colorTrack); 
   x = tracksCoordAmpLink(i,1:8:end); 
   y = tracksCoordAmpLink(i,2:8:end); 
   clear color idx D
   
%    if iFrame == numel(analInfo);
%        color = 'k';
%    else
%        currentFrame = 4+(iFrame-1)*8;
%        nextFrame = 4+(iFrame)*8;
%        % figure out the local frame to frame velocity
%        vel = diff(tracksCoordAmpLink(i,[currentFrame,nextFrame]));
%        if isnan(vel)
%            color = 'k';
%        else
%            vel = vel*(pixSize_um/timeInt);
%            
%            D=createDistanceMatrix(vel,mapper);
%            [sD,idx]=sort(abs(D),2);
%            idx=idx(:,1);
%            
%            color = cell2mat(cMap(idx,:));
%            % plot the current frame bigger
 %       end 
% color = 'k';
%       scatter(x(iFrame),y(iFrame),50,color,'filled');
       
 %  end % if frame ==1
   end % if there is a track
  clear color idx D vel
   text(xCoordsPostFilt(i,iFrame),yCoordsPostFilt(i,iFrame),num2str(trackID(i)),'color','k'); 
   text(10,ny-10,num2str(iFrame),'color','k');
end 

  cellfun(@(x) plot(x(:,2),x(:,1),'k'),roiYXNE); 

 print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
     [saveDir filesep 'frame_' num2str(iFrame, fmt) ext])
 saveas(h,[saveDir filesep 'frame_' num2str(iFrame,fmt) '.fig']); 
 % save a frame 
 if iFrame == 101
    saveas(h,[saveDir filesep 'frame_' num2str(iFrame,fmt) '.eps'],'psc2'); 
  end
  close(gcf)
end 
%% Make avi Movie
cd(saveDir)
execute = 'mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
system(execute); 
%% plot time series 



   

%%
lengthPlots = 1; 
if lengthPlots == 1
 if ~isdir([saveDir filesep 'timeSeries' filesep 'lengths']) 
        mkdir([saveDir filesep 'timeSeries' filesep 'lengths'])
        mkdir([saveDir filesep 'timeSeries' filesep 'velocities'])
 end 
for iTrack = 1:length(trackID)
    % extract amps (ie lengths)
    lengths = tracksCoordAmpLink(iTrack,4:8:end); 
   % lengths = lengths(~isnan(lengths));
if lengthPlots ==1 
  % reset handle
 %   h = figure('Visible','off','Position',[50 50 nx ny]) ;

%iptsetpref('ImshowBorder','tight');
% 
% set(h, 'InvertHardcopy', 'off');
% 
% set(h, 'PaperUnits', 'Points');
% 
% set(h, 'PaperSize', [nx ny]);
% 
% set(h, 'PaperPosition', [0 0 nx ny]); % very important
% 
%  set(h,'DefaultLineLineSmoothing','on');
%  set(h,'DefaultPatchLineSmoothing','on');
% 
% Configure axes 
h = figure('Visible','on'); 
fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];

    
% Create axes
axes1 = axes('Parent',h,'FontWeight','bold','FontSize',14,...
    'FontName','Arial');
hold(axes1,'all');
%lengths = lengths.*pixSize_um; 
y = 0:length(lengths)-1; 
    y = 5.*y/60; 
% Create plot
plot(y,lengths,'Parent',axes1,'MarkerFaceColor','g','Marker','o','MarkerSize',5,'Color','b');

% Create xlabel
xlabel({'Time (Min)'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');

% Create ylabel
%ylabel({'Length','um'},'FontWeight','bold','FontSize',14,...
 %   'FontName','Arial');
 ylabel({'Veil Velocity','nm/sec'},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');

% Create title
%title({'Time Series: Filopodia Length (um)';['Track' num2str(trackID(iTrack))]},'FontWeight','bold','FontSize',14,...
    %'FontName','Arial');
  
title({'Time Series: Veil Velocity (nm/sec)';['Track' num2str(trackID(iTrack))]},'FontWeight','bold','FontSize',14,...
    'FontName','Arial');

%      print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
%      [saveDir filesep 'timeSeries' filesep 'lengths' filesep 'track_' num2str(iTrack, fmt) ext])

%saveas(h,[saveDir filesep 'timeSeries' filesep 'lengths' filesep 'track_' num2str(trackID(iTrack), fmt) '.eps'],'psc2'); 
%saveas(h,[saveDir filesep 'timeSeries' filesep 'lengths' filesep 'track_' num2str(trackID(iTrack), fmt) '.fig'])
saveas(h,[saveDir filesep 'timeSeries' filesep 'lengths' filesep 'track_' num2str(trackID(iTrack), fmt) '.fig'])
close(gcf) 
end % if length
%% 
velPlots = 1; 
if velPlots ==1
    % PLOT VELOCITY TIME SERIES %
    % reset handle
    h = figure('Visible','off','Position',[50 50 nx ny]) ;
    
    %iptsetpref('ImshowBorder','tight');
    
    set(h, 'InvertHardcopy', 'off');
    
    set(h, 'PaperUnits', 'Points');
    
    set(h, 'PaperSize', [nx ny]);
    
    set(h, 'PaperPosition', [0 0 nx ny]); % very important
    
    % set(h,'DefaultLineLineSmoothing','on');
    % set(h,'DefaultPatchLineSmoothing','on');
    
    % Configure axes
    
    fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];
    % Create axes
    axes1 = axes('Parent',h,'FontWeight','bold','FontSize',14,...
        'FontName','Arial');
    hold(axes1,'all');
    y = 0:length(lengths)-2; 
    y = 5.*y/60; 
    
    % Create plot
    plot(y,diff(lengths,1,2).*(pixSize_um/timeInt*60),'Parent',axes1,'MarkerFaceColor','r','Marker','o','MarkerSize',5,'Color','b');
    
    % Create xlabel
    xlabel({'Time (Min)'},'FontWeight','bold','FontSize',14,...
        'FontName','Arial');
    
    % Create ylabel
    ylabel({'Velocity','um/min'},'FontWeight','bold','FontSize',14,...
        'FontName','Arial');
    
    % Create title
    title({'Time Series: Filopodia Velocities ';['Track' num2str(trackID(iTrack))]},'FontWeight','bold','FontSize',14,...
        'FontName','Arial');
    
    
    %       print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
    %      [saveDir filesep 'timeSeries' filesep 'velocities' filesep 'track_' num2str(iTrack, fmt) ext])
   % saveas(h,[saveDir filesep 'timeSeries' filesep 'velocities' filesep 'track_' num2str(trackID(iTrack), fmt) '.eps'],'psc2')
    saveas(h,[saveDir filesep 'timeSeries' filesep 'velocities' filesep 'track_' num2str(trackID(iTrack), fmt) '.fig'])
   close(gcf)
    
end % if plotVel



end % iTrack


outputTracking = plotHistsVelocityFilo(tracksCoordAmpLink,saveDir); 


end 