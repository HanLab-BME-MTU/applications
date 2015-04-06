function [neuriteLength,error,pathMasksFinal] = GCAfindNeuriteLength( analInfo ,saveDir,imDir,plotPathColor,paramsIn)
% GCAfindNeuriteLength: reads in the the output of GCA analyzer detection
% and calculates the global parameter neurite length in each frame. 
% The neurite length is defined here as the longest path from the neurite 
% entrance point into the frame to the tip of the neurties veil/stem skeleton 
%(Note to self : might
% want to make this a single frame input instep of doing the wrapper per
% frames here... ) you also might want to take out the plotting functions
% and do that after make sure the output is ameniable to plotting as well)
% originally getNeuritePathForFinal (renamed 20140721)
%
% INPUT:
% analInfo: the output of the neurite segmentation for a given movie (NOTE may in the future
%           separate the output so that reads in the body estimate separately so
%           saves space
% mkPlot: (scalar) 1 if mkPlot option on
% saveDir: (char) path For Saving
% imDir: (char) path For Images % this will later be movieData formate
% plotting structure 
% have option for verbose output versus streamlined output 
% verbose output will include the 
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% START INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % INPUT R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% analInfo: will shift this around but basically you need the bodyMask
% wherever that is stored and some pointers to the image directory and that
% is it. 
% 
% saveDir: where are you going to store your output
% 
% imDir : again this can be eventually taken care of in movieData 
% 
% 
% 
%     params are going to be put into 1 structure called paramsIn           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INPUT: PARAMS: MOVIE OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% makePresentMovie: 1 will make an output directory 
%                  'Present/PresentOverlays/' containing the list of .png
%                  files for the simple movie with overlay 
% 
% plotAllPaths: 1 will make a an outputdirectory 
%               'TroubleshootMovieFile'
%               'All Paths'
%
% plotMoviePathsWithDist: 1 wll make an outputdirectory 
% 'TroubleshootMovieFile' 
%    'All Paths Plus DistTrans'
%
% plotPathColor : 'default will be 'g'
%
% neuriteBodyColor: default will be 'y'
%
% makeScatterMovie: 1 will make an output directory '
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       PARAMS: VERBOSITY OPTIONS  
%(ie how much of the information regarding the movie you would like to 
% calculate and store)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% savePathTifs: 1 will make an output directory 
%   LongPathTifs : with a logical overlay of the longest path .tif files 
%   recommeded if you ever need to replot anything without rerunning the
%   algorithm Default = 1 
%   
% 
% firstLastFrameOnly: if 1 will save time to just get a net outgrowth
%               parameter useful for screening neurites quickly or for neurites that
%               might have semgmentation problems mid movie but still have useful data thus
%               need an outgrowth parameter 
%               Default = 0 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Paths : file of .tifs marking the paths again these will be called by teh
% the movie making functions so the user can make a movie at a later date
% if they save this file without re-running the original function  - set up
% currently as .tif files (could maybe store as a structure of 
%  .xyCoords or the pixIndices
% in order - you would have to convert them back to  
% 
% 
%  
% 
%% Check/Set Parameters  

if nargin < 4 % switch after
    paramsIn.plotPathColor ='b' ; %PUT INTO PARAMETERS! 
else 
    paramsIn.plotPathColor = plotPathColor; % quick fix 20140216- now can 
    % the usefuness of the 'params' as you can quickly change 1 think and
    % have the rest default. fix your other functions before release. ... 
end 

if nargin<5 
    
    
    paramsIn.OutputDirectory = [saveDir  filesep 'PARAMETER_EXTRACTION' filesep 'GlobalFunctional' filesep  'neurite_outgrowth_measurements'];
    paramsIn.OutputDirectoryBodyShape = [saveDir filesep 'PARAMETER_EXTRACTION' filesep 'Veil' filesep 'veil_size_measurements']; % for now I know this is 
    % a bit weird to  
    paramsIn.makePresentMovie = 1;
    paramsIn.makeScatterMovie = 0;
    paramsIn.plotAllPaths = 1;
    paramsIn.plotMoviePathsWithDist = 1;
    %paramsIn.plotPathColor ='m' ;
    paramsIn.neuriteBodyColor =  'y';
    paramsIn.longPathTifs = 1;
    paramsIn.firstLastFrameOnly = 0;
    paramsIn.makeMovieFiles = 0;
    paramsIn.distCutOffForBody = 30;  % currently in um
    paramsIn.mkPlotsBodyMeasure = 1;
    paramsIn.subRegions =1; % flag to make frontMask subRegions
    paramsIn.plotFrontMask = 1; % flag to plot the frontMask 
end
p = paramsIn; 
%% Make all output files 
mkClrDir(p.OutputDirectory)
% if make movies initiate a cell array for the movie directories 
if p.makeMovieFiles ==1 
    numMovies = (p.makePresentMovie + p.makeScatterMovie + p.plotAllPaths + p.plotMoviePathsWithDist) ; 
    movieDirs = cell(numMovies,1);
end 
if p.makePresentMovie  == 1 
  presentDirOverlay =  [p.OutputDirectory filesep 'Present' filesep 'PresentOverlays']; 
  mkClrDir(presentDirOverlay); 
 
  if p.makeMovieFiles == 1
  emptyDirs  = cellfun(@(x) isempty(x), movieDirs); 
  toFill = find(emptyDirs,1,'first');  
  movieDirs{toFill} = presentDirOverlay; 
  end 
end 

if p.makeScatterMovie  == 1
    scatterDir = [p.OutputDirectory filesep 'Present' filesep 'ScatterMovie']; 
    mkClrDir(scatterDir)
    if p.makeMovieFiles == 1
        emptyDirs  = cellfun(@(x) isempty(x), movieDirs);
        toFill = find(emptyDirs,1,'first');
        movieDirs{toFill} = scatterDir;
    end
end 

if p.plotAllPaths == 1
    allPathDir = [p.OutputDirectory filesep 'Troubleshoot' filesep 'AllPaths'];
    mkClrDir(allPathDir)
    if p.makeMovieFiles == 1
        emptyDirs  = cellfun(@(x) isempty(x), movieDirs);
        toFill = find(emptyDirs,1,'first');
        movieDirs{toFill} = allPathDir;
    end
end 

if p.plotMoviePathsWithDist == 1
    pathsWithDist = [p.OutputDirectory filesep 'Troubleshoot' filesep 'AllPathsWithDist']; 
    mkClrDir(pathsWithDist); 
    if p.makeMovieFiles == 1
        emptyDirs  = cellfun(@(x) isempty(x), movieDirs);
        toFill = find(emptyDirs,1,'first');
        movieDirs{toFill} = pathsWithDist;
    end
end 

if p.longPathTifs == 1 
    pathTifDir = [p.OutputDirectory filesep 'Output' filesep 'LongPathMasks']; 
    mkClrDir(pathTifDir); 
end 

if p.subRegions == 1
    subRegDir = [p.OutputDirectory filesep 'Output' filesep 'frontMasks']; 
    mkClrDir(subRegDir); % need to make a directory to save the output masks 
    subNames{1} = 'Front'; 
    subNames{2} = 'Back'; 
    subNames{3} = 'Side'; 
    
      arrayfun(@(i) mkClrDir([subRegDir filesep subNames{i} filesep 'masks']),1:3);   
    if p.plotFrontMask == 1
            plotFrontMaskDir = [p.OutputDirectory filesep 'Output' filesep 'frontMasks' filesep 'Troubleshoot']; 
            mkClrDir(plotFrontMaskDir); 
    end 
end 

%% 
error = 0;
neuriteLength = zeros(numel(analInfo),1);
if paramsIn.firstLastFrameOnly == 1
    frames = [1,numel(analInfo)];
else
    frames = 1:numel(analInfo);
end

for i = 1:length(frames)
    iFrame = frames(i);
    display(['Calculating Longest Distance Path for Frame' num2str(iFrame)]); 
    mask = analInfo(iFrame).masks.neuriteEdge;
    [ny,nx] = size(mask); 
    roiYX  = bwboundaries(mask);
%% peform thinning     
   % thinnedBodyAll = bwmorph(mask,'thin','inf');
   % I think that the 2D skeletonization method is more appropriate here. 
   thinnedBodyAll = bwmorph(mask,'skel','inf'); 

   
        
  
%% CREATE GRAPH FROM SKELETON 
   [adjMat,edges , tipVertices,tipXYCoords,vertMat,edgePixels]   = skel2Graph4OutGrowthMetric(thinnedBodyAll);
   
     %% Small fix added 2014-05-04 to show different path examples 
 forPresent = 1;
  if forPresent ==1
     
      listOfImages = searchFiles('.tif',[],imDir,0);
        img = double(imread([listOfImages{iFrame,2} filesep listOfImages{iFrame,1}]));
        [ny,nx ] = size(img); 
       setFigure(nx,ny,'off'); 
        
        imshow(-img,[])
       hold on 
       spy(thinnedBodyAll); 
cmap = hsv(numel(edgePixels)); 
       
   for iEdge = 1:numel(edgePixels) 
     
      
       
    
[y,x]  = ind2sub(size(img),edgePixels{iEdge});

 plot(x,y,'color',cmap(iEdge,:), 'Linewidth',2); 

 
   end 
      nVerts = length(unique(edges(:))); 
   for iVert = 1:nVerts
        [y,x] = ind2sub(size(thinnedBodyAll),find(vertMat==iVert));
        text(x(1),y(1),num2str(iVert),'color','y');
        hold on
    end
   
   cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX); 
    saveas(gcf,['testColoredPaths' num2str(iFrame,'%03d') '.eps'],'psc2')
 close gcf
  end 
 %%   
   
   
   
   
   
   
   
   
%    fix =0 ; 
%    if fix == 1
%        % load backbone info 
%     load([saveDir filesep 'backboneInfo.mat']);
%    [coords] = backboneInfo(iFrame).coordsEnterNeurite; 
%    xEnter  = coords(1,1); 
%    yEnter = coords(1,2); 
   
%    else  (Not sure what the point of this little code here was- can
%    likely remove as I was automatically setting fix =0 
   
   %%% FIND THE SKELETON TIP THAT IS CLOSEST TO THE NEURITE COORD %%%%  
   idxEnter  = analInfo(iFrame).idxEnterNeurite;
   
   [yEnter,xEnter] = ind2sub(size(thinnedBodyAll),idxEnter);
%    end 
   nTips = length(tipXYCoords(:,1));
   
   %get distances of each tip to the enter coord 
   distFromEnt = arrayfun(@(i) sqrt((xEnter-tipXYCoords(i,1))^2 + (yEnter-tipXYCoords(i,2))^2),1:nTips);
   
   % find tip with teh min distance: this will be start point for the graph
   % shortest path algorithm 
   enterVertex = tipVertices(distFromEnt==min(distFromEnt));
   
   % the rest of the skeleton tips will be the other input into the graphshortestpath algorithm 
   tipVertices(distFromEnt==min(distFromEnt)) = []; % take out the enter index
   
   %Check if tipVertices empty: this indicates something is wrong.
   %which shouldn't happen in a perfect world but who are we kidding..
   %if there is problem flag it and just use the paths from the previous
   %frame for all the calcs and the overlays so things don't completely
   %crash for the movie
   if isempty(tipVertices)
       display(['Morphology of skeleton weird for Frame ' num2str(iFrame) ': using Frame' num2str(iFrame-1)])
       iFrame = iFrame-1;
       error = 1;
       mask = analInfo(iFrame).masks.neuriteEdge;
       [ny,nx] = size(mask); 
       roiYX  = bwboundaries(mask);
       
       thinnedBodyAll = bwmorph(mask,'thin','inf');
       [adjMat,edges , tipVertices,tipXYCoords,vertMat,edgePixels]   = skel2Graph4OutGrowthMetric(thinnedBodyAll,0);
       
       %   plot the shortest path from tipVertices to the idxEnterNeurite
       
       % find the path from the tip end point
       
       idxEnter  = analInfo(iFrame).idxEnterNeurite;
       
       [yEnter,xEnter] = ind2sub(size(thinnedBodyAll),idxEnter);
       
       % find the tip that is closest to the pointof entry
       nTips = length(tipXYCoords(:,1));
       distFromEnt = arrayfun(@(i) sqrt((xEnter-tipXYCoords(i,1))^2 + (yEnter-tipXYCoords(i,2))^2),1:nTips);
       enterVertex = tipVertices(distFromEnt==min(distFromEnt));
       tipVertices(distFromEnt==min(distFromEnt)) = []; % take out the enter index 
   end % isempty(tipVertices)
   
%% FIND THE LONGEST PATH OF THE SKELETON FROM THE NEURITE ENTRANCE TO THE END OF THE GROWTH CONE
    % This function will find the shortest path between the neurite enter vertex
    % and each of the tips of the neurite skeleton on the growth cone given the
    % distance edge weights 
    % from there one can choose the longest of these paths as a good approximation of the neurite 
    % length metric for a given frame 
    for iTip = 1:length(tipVertices)
        [distPath(iTip), path{iTip} ] = graphshortestpath(adjMat,enterVertex,tipVertices(iTip),'Directed',false);
    end
    
    if forPresent == 1 
        imshow(-img,[])
        
    end 
        
   
        % (MY NOTES TO EDIT)  If have singletons I think sometimes the distPath goes to inf- before
        % when I had the length metric not converted to xy coords (in raw
        % pixels) I must have had some filtering metric for singletons that 
        % was eliminated when I fixed the length metric --- check. 
        % 2014_01_04 
        % take out empty paths here and see if it is a quick fix 
        idxOut = cellfun(@(x) ~isempty(x),path); 
        path =  path(idxOut); 
        distPath = distPath(idxOut); 
        
         if isempty(distPath)
     
        % something is wrong MARK! 
        display(['Could NOT find the longest Path for frame ' num2str(iFrame,'%03d')])
        display(' Filling Neurite Global Path With NaN'); 
        neuriteLength(iFrame) = NaN; 
    else 
        

        % Find max distance path 
        idxMax =  find(distPath==max(distPath));
        idxMax  = idxMax(1); % just choose the first in case there is a tie.
        longPath = path{idxMax};
%% MAKE FINAL LONG PATH MASKS AND RECALCULATE THE FINAL DISTANCE OF THE LONG PATH 
% NOTE: in the original implementation the length to the vertices is not included in the final 
%       distPath calc therefore it is slightly more correct to remake the path mask
%       with the vertices and re-calculate the distance
%       -20140425 should be fixed check this maria and give final ok before
%       release 

        edgesAll = [vertcat(edges(:,1), edges(:,2)), vertcat(edges(:,2),edges(:,1))];
        idxEdge = zeros(length(longPath)-1,1);
        % plot the pixels that correspond to the path lenghs
        for iPath = 2:length(longPath)
            startNode = longPath(iPath-1);
            endNode = longPath(iPath);
            %test(i-1,:) = [startNode,endNode];
            idxLog = arrayfun(@(x) isequal(edgesAll(x,:),[startNode,endNode]),1:length(edgesAll(:,1)));
            idxEdge(iPath-1) = find(idxLog==1);
        end
        edgePixels = [edgePixels';edgePixels'] ;
        
        
        % get edge Pixels
        pathMask = zeros(size(thinnedBodyAll));
        pathMask(vertcat(edgePixels{idxEdge})) = 1;
        

         idxVert = arrayfun(@(x) find(vertMat==x),longPath,'uniformoutput',0); 
          % add the path vertices to the path mask 

         pathMask(vertcat(idxVert{:})) = 1; 
       
        thinMask = 1; 
while thinMask == 1
    pathMask   = bwmorph(pathMask,'thin','inf'); 
    pathMask = bwmorph(pathMask,'spur'); 
     pathMask(vertcat(idxVert{[1,end]})) = 1; % make sure to add back the vertices 
    
    % keep testing the path mask for junctions 
        nn = padarrayXT(double(pathMask~=0), [1 1]);
        sumKernel = [1 1 1];
        nn = conv2(sumKernel, sumKernel', nn, 'valid');
        nn1 = (nn-1) .* (pathMask~=0);
     
        test= nn1 > 2;
      % no more problem points you are done 
        if sum(test(:)) ==0 
         thinMask = 0; 
        end 
end 
        
%% added 20140426 %% 



    if p.longPathTifs ==1 
        imwrite(pathMask, [pathTifDir filesep num2str(iFrame,'%03d') '.tif']);  
    end 
         
         
         
         
         finalPathPix= find(pathMask==1);
         EPFinalPath = getEndpoints(finalPathPix,[ny,nx],1); 
         [nEndpoints,~]= size(EPFinalPath); 
         dist = arrayfun(@(i) sqrt((xEnter-EPFinalPath(i,1))^2 + (yEnter-EPFinalPath(i,2))^2),1:nEndpoints); 
         
         
         EPFinalPath = EPFinalPath(dist==max(dist),:); % take the end farthest away from the entrance
         
       % Reorder the pixels 
        pixIdxBack = nan(length(finalPathPix),1); % overinitialize to make happy     
        transform = bwdistgeodesic(pathMask,EPFinalPath(:,1),EPFinalPath(:,2)); % input in col(x),row(y) (xy coords)
        iPix = 0; 
           while length(find(transform==iPix)) == 1
              pixIdxBack(iPix+1) = find(transform==iPix); % start at the endpoint
              iPix = iPix +1;  
           end   
        pixIdxBack = pixIdxBack(~isnan(pixIdxBack)); 

        % Recalculate the distance 
       [ neuriteLength(iFrame),pixForBodyEst]  = calculateDistance(pixIdxBack,[ny,nx],0,'distCutoff',p.distCutOffForBody);% output in um
       % perform distance transform from the VEIL/STEM border 
       
       if p.subRegions ==1; 
       distTrans = bwdist(~mask); 
       distVals = distTrans(pixForBodyEst); 
       subRoiInfo = [pixForBodyEst distVals];
      [subRois,xVect,yVect] =   GCAgetGrowthConeSubRegions(subRoiInfo,45,[ny,nx],4,mask); 
     % save in a mask file for later use: Think about best place for
     % substructure to remain in line with movie data 
     
     arrayfun(@(i) imwrite(subRois(:,:,i),[subRegDir filesep subNames{i} filesep 'masks' filesep 'mask' num2str(iFrame,'%03d') '.tif']),1:3); 
     
       
       end  
      
       
       
       end 
          % if 
          
          
            % to fill these in it might be just easiest to use bresenham here
            % as otherwise need to get the vertices index again... 
            pathMasksFinal(:,:,iFrame) = pathMask;

%% Optional Plots : troubleshoot 
        rotate =0; 
if (p.plotAllPaths ==1 || p.makePresentMovie)
   
        listOfImages = searchFiles('.tif',[],imDir,0);
        img = double(imread([listOfImages{iFrame,2} filesep listOfImages{iFrame,1}]));
     if rotate == 1
         angle = 90;
        img = imrotate(img,angle); 
        
     end 
     
     %   imgLarge = [img img];
     
            
      %  [nyLarge,nxLarge] = size(imgLarge);
        [ny,nx] = size(img); 
      %  h = (nxLarge,nyLarge,'off');
h = setFigure(nx,ny,'off') ; 
        imshow(-img,[]) ;
        hold on
        if rotate == 1
           mask  = imrotate(mask,angle); 
           roiYX =  bwboundaries(mask); 
           pathMask = imrotate(pathMask,angle);
           finalPathPix = find(pathMask ==1); 
           EPFinalPath = getEndpoints(finalPathPix,[ny,nx],1); 
         EPFinalPath = EPFinalPath(1,:); % just take the first for ordering  
           % Reorder the pixels 
        pixIdxBack = nan(length(finalPathPix),1); % overinitialize to make happy     
        transform = bwdistgeodesic(pathMask,EPFinalPath(:,1),EPFinalPath(:,2)); % input in col(x),row(y) (xy coords)
        iPix = 0; 
           while length(find(transform==iPix)) == 1
              pixIdxBack(iPix+1) = find(transform==iPix); % start at the endpoint
              iPix = iPix +1;  
           end   
        pixIdxBack = pixIdxBack(~isnan(pixIdxBack)); 

        end 
        cellfun(@(x) plot(x(:,2),x(:,1),'color','y','Linewidth',2),roiYX);
        hold on
        % plot the long path
        % spy(pathMask,p.plotPathColor,10)
        [yPath,xPath]  =  ind2sub([ny,nx],pixIdxBack);
        plot(xPath,yPath,'color',p.plotPathColor,'Linewidth',2); 
   
        % Add the bells and whistles 
    pixSizeMic = 0.216; % in um
    timeInterval = 5; % in sec
    
   % text(5,ny-20,[{'Neurite Length = ';[ '         ' num2str(neuriteLength(iFrame),3) ' um']}],'color','k','FontSize',12,'FontName','Arial','FontWeight','bold');
     %text(5,ny-30,'','Color','k','FontSize',12,'FontName','Arial','FontWeight','Bold'); 
    %text(5,ny-10,[ num2str((iFrame-1).*timeInterval) ' sec'],'color','k'); % add variable of time interval
    %pixels = round(10/pixSizeMic);  
    %plotScaleBar(pixels,pixels/10,'Color',[0 0 0],'location','southeast');
    
    if p.makePresentMovie == 1
        saveas(gcf,[presentDirOverlay filesep num2str(iFrame,'%03d') '.png']); 
        saveas(gcf,[presentDirOverlay filesep num2str(iFrame,'%03d') '.eps'],'psc2');
    end 
    cd(presentDirOverlay)
   
    % add rest 
    if p.plotAllPaths==1 
    close gcf % for now just make a new figure 
    setFigure(nx,ny,'off');
    imshow(-img,[]);
    hold on 
    spy(thinnedBodyAll,'b')
    cellfun(@(x) plot(x(:,2),x(:,1),'y'),roiYX); 
    hold on
    %spy(pathMask,p.plotPathColor,10); % see if can erase from previous figure and then reoverlay 
    plot(xPath,yPath,'color',p.plotPathColor,'Linewidth',2);
    nVerts = length(unique(edges(:))); 
    % plot the vertice with it's label
    for iVert = 1:nVerts
        [y,x] = ind2sub(size(thinnedBodyAll),find(vertMat==iVert));
        text(x(1),y(1),num2str(iVert),'color','y');
        hold on
    end
    
    text(5,ny-30,['Neurite Length = ' num2str(neuriteLength(iFrame),3) 'um'],'color','k');
    text(5,ny-10,[ num2str((iFrame-1).*timeInterval) ' sec'],'color','k'); % add variable of time interval
    pixels = round(10/pixSizeMic);  
    plotScaleBar(pixels,pixels/10,'Color',[0 0 0],'location','southeast');
    
    
    saveas(gcf,[allPathDir filesep num2str(iFrame,'%03d') '.png']); 
    saveas(gcf,[allPathDir filesep num2str(iFrame,'%03d'),'.fig']); 
    
    
    end 
      
    close gcf
end % if (p.plotAllPaths == 1 || p.makePresentMovie)
    
if p.plotMoviePathsWithDist == 1
    % plot the dist trans
    distTrans = bwdist(~mask);
    distTrans = distTrans.*pixSizeMic;
    setFigure(nx,ny,'off');
    overlay = distTrans.*thinnedBodyAll;
    imagesc(overlay);
    hold on
    cellfun(@(x) plot(x(:,2),x(:,1),'color','w'),roiYX);
    colorbar
    saveas(gcf,[pathsWithDist filesep num2str(iFrame,'%03d') '.png']);
    saveas(gcf,[pathsWithDist filesep num2str(iFrame,'%03d') '.fig']); 
    close gcf
end

if ~isempty(p.distCutOffForBody); 
   if ~isdir(p.OutputDirectoryBodyShape)
       mkdir(p.OutputDirectoryBodyShape); 
   end 
   cd(p.OutputDirectoryBodyShape);
   
    % save the body shape values 
    distTrans = bwdist(~mask); 
    distTrans = distTrans.*pixSizeMic; 
   
    
    valuesBody= distTrans(pixForBodyEst); 
    findMaxDist  = 0; 
 % Begin to get the thickest part of the GC for potential regional metrics
 if findMaxDist == 1
   [yCenter,xCenter] =  ind2sub([ny,nx],pixForBodyEst(valuesBody==max(valuesBody))); 
 else
     [yCenter,xCenter] = ind2sub([ny,nx],pixForBodyEst(end)); 
     
 end 
    bodyShape{iFrame} = valuesBody; 
     if p.mkPlotsBodyMeasure == 1 
        visBodyEstDir = [p.OutputDirectoryBodyShape filesep 'Visuals']; 
        if ~isdir(visBodyEstDir)
            mkdir(visBodyEstDir); 
        end
        
        fsFigure(.75,'visible','off');
        
        subplot(1,2,2);
       tipGCMask = zeros([ny,nx]); % this the path from the tip of the longest path to  a user selected distance along the growth cone
       tipGCMask(pixForBodyEst) = 1; 
       % get pixels for body estimation from the medial axis transform 
    
       [yForBody,xForBody] = ind2sub([ny,nx],pixForBodyEst); 
       imagesc(tipGCMask.*distTrans,[0 5])
    
       hold on 
       roiYX = bwboundaries(mask); 
       cellfun(@(x) plot(x(:,2),x(:,1),'w'),roiYX); 
       set(gca,'xTickLabel',[]); 
       set(gca,'yTickLabel',[]); 
       pixels = round(5/pixSizeMic);  
       plotScaleBar(pixels,pixels/10,'Label','5 um','Color',[1 1 1]);
          scatter(xCenter,yCenter,100,'w','filled')
       colorbar
      % title({'AvgDistToBoundary' num2str(mean(bodyShape{iFrame}),3) 'um' ; 'Along'  num2str(p.distCutOffForBody) ' um From Tip'})
       subplot(1,2,1)
       
       imagesc(img)
       hold on 
        plot(xForBody,yForBody,'w','Linewidth',2);
        colorbar 
      
       cellfun(@(x) plot(x(:,2),x(:,1),'w'),roiYX);  
        title('Image')
       set(gca,'xTickLabel',[]); 
       set(gca,'yTickLabel',[]); 
       
       
       saveas(gcf,[visBodyEstDir filesep num2str(iFrame,'%03d') '.fig']); 
       saveas(gcf,[visBodyEstDir filesep num2str(iFrame,'%03d') '.png']); 
       
       close gcf
     end 
     
if (p.subRegions && p.plotFrontMask) 
    c = ['r','g','b'];  % for now front is red, green is back, and blue is sides. 
    
    setFigure(nx,ny,'off'); 
       imshow(-img,[]) ; 
       hold on 
       % plot the subrois 
       for iSub = 1:3 
           
       roiYXSub = bwboundaries(subRois(:,:,iSub)); 
       
       cellfun(@(x) plot(x(:,2),x(:,1),'color',c(iSub),'Linewidth',3),roiYXSub); 
       end % iSub 
       [yFront,xFront] = ind2sub([ny,nx],pixForBodyEst); 
       
       %cmap = jet(length(xFront)); 
       %arrayfun(@(i) scatter(xFront(i),yFront(i),10,cmap(i,:,:),'fill'),1:length(xFront)); 
      % color by distTrans 
%        cMapLength = 128; cMap= jet(cMapLength); % create the map- could maybe make so have more colors
%             distVals = subRoiInfo(:,2); 
%             % convert to um 
%             distVals = distVals.*0.216;
%             % set upper limits to 10 um
%             distVals(distVals>10) = 10; 
            % mapper
%             mapper=linspace(0,10,cMapLength)'; % CURRENTLY PLOTTING FROM 0 to 20 um/min
          
             % get closest colormap index for each feature
%              D=createDistanceMatrix(distVals,mapper);
%              [sD,idx]=sort(abs(D),2);
%              
%              
%              % plot according to dist trans 
%              arrayfun(@(k) scatter(xFront(idx(:,1)== k),yFront(idx(:,1)==k), 10, cMap(k,:),'filled'),1:cMapLength);  
             
%               x=xMatIn(:,1:end-1)';% make so row is time frame, col is subtrack number 
%               y=yMatIn(:,1:end-1)';
%              for k=1:cMapLength
%                 plot(x(:,idx(:,1)==k),y(:,idx(:,1)==k),'color',cMap(k,:),'lineWidth',1);
%              end
       
       
        u = xVect(1) - xVect(end); % vector towards the tip of neurite
        v = yVect(1) - yVect(end); 
       
        % Plot the direction of the local vector in black 
        
        quiver(xVect(end),yVect(end),u,v,10,'color','k','Linewidth',2); 
        
        plot(xFront,yFront,'color','k','linewidth',1); 
        % scatter the points used in the calculation of local vector. 
        scatter(xVect(1),yVect(1),50,'k','fill');
        scatter(xVect(end),yVect(end),50,'k','fill'); 
        % scatter the center point around which the local vector was
        % plotted
        scatter(xCenter,yCenter,100,'w','fill','p'); 
        
        % give a title 
        
        saveas(gcf,[plotFrontMaskDir filesep 'frontMaskPlot' num2str(iFrame,'%03d') '.fig']);
        saveas(gcf,[plotFrontMaskDir filesep 'frontMaskPlot' num2str(iFrame,'%03d') '.png']);
       close gcf
    
end 
     
    
 


% if (iFrame == 1 && p.pathMovie ==1 ) 
%     c = hsv(numel(paths)); 
%     for iPath = 1:numel(paths)
%     setFigure(nx,ny,'on')
%     
%     imshow(-img,[]);
%     hold on 
%     
%     spy(thinnedBodyAll,'b');
%     path{iPath} 
%     
    
  %  end 
    
    
    
    clear longPath distPath path pathMask
    
    end % if something fucked up 
end % iFrame = 1:lengthFrames 

%% SAVE OUTPUT 
timeStamp = clock; % for now keep as a separate variable eventually link! 
save([p.OutputDirectory filesep 'neuriteLengthOutput.mat'],'neuriteLength','timeStamp');


% save body calcs if applicable 
if ~isempty(p.distCutOffForBody)
    
      save([p.OutputDirectoryBodyShape filesep 'bodyShape.mat'],'bodyShape','timeStamp'); 
  
end 
 
%% MAKE SCATTER MOVIE (Optional) 
if p.makeScatterMovie
    plotNeuriteOutgrowthInTime(neuriteLength,p.plotPathColor,1,timeInterval,0,scatterDir,'single')   
end 
%% FINISH MAKING MOVIE FILES (OPTIONAL)   
if p.makeMovieFiles == 1
    for iMovie = 1:numel(movieDirs)
        cd(movieDirs{iMovie})
        execute = 'mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
        system(execute);
    end
end
end

