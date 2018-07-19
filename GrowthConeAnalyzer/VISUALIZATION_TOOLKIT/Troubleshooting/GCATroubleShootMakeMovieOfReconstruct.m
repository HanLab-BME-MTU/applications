function [ hSet,filoFilterSet,filoParams ] = GCATroubleShootMakeMovieOfReconstruct(filoBranch,veilStem,frame,pixSizeMic,imDir,varargin)
% Small function that makes a nice movie of the main reconstruction steps
% As the name suggests it is useful for troubleshooting segmentation errors
%% check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addParameter('writeTitles',true);

ip.addParameter('plotScaleBar',true);
ip.addParameter('sizeScaleBar', 10); % in um 

ip.parse(varargin{:});
%%
listOfImages = searchFiles('.tif',[],imDir);
filename = [char(listOfImages(frame,2)) filesep char(listOfImages(frame,1))];
textColor = [ 0 0 0 ];
fontText =  {'FontName','Arial','FontSize',14,'FontName','Arial','color',textColor};
%%%% START %%%%
img = double(imread(filename));

[ny,nx] = size(img);

if ip.Results.plotScaleBar
     pixels = 10/pixSizeMic;
end 

%% Initiate

countFig = 1;
%% 00 -original image

hSet(countFig).h =setFigure(nx,ny);
imshow(-img,[]) ;
hold on

if ip.Results.writeTitles
    text(nx/10, 10,'Original Image ', fontText{:});
end

if ip.Results.plotScaleBar
    plotScaleBar(pixels,pixels/10,'Color',textColor);
end
countFig = countFig+1;

%% Final Body Mask %%%%
hSet(countFig).h = setFigure(nx,ny);
imshow(-img,[]) ;
hold on
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);

cellfun(@(x) plot(x(:,2),x(:,1), 'color',[0.0039  ,  0.4264 ,   0.3848],'Linewidth',1),edgeYX);
if ip.Results.writeTitles
    text(nx/10, 10,'Veil/Stem Estimation Complete', fontText{:});
end

if ip.Results.plotScaleBar
    plotScaleBar(pixels,pixels/10,'Color',textColor);
end
countFig = countFig+1;
%% Show Seed
hSet(countFig).h = setFigure(nx,ny);
imshow(-img,[]) ;
hold on

if ip.Results.writeTitles
    text(nx/10, 10,'Get Seed For Reconstruction', fontText{:});
end

seedMask = filoBranch(frame).reconstructInfo.seedMask{1};
spy(seedMask,'b');

if ip.Results.plotScaleBar
    plotScaleBar(pixels,pixels/10,'Color',textColor);
end

countFig = countFig+1;
%%  Show Candidates
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
if ip.Results.writeTitles
    text(nx/10, 10,'Get Candidates', fontText{:});
end
spy(seedMask,'b');
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
preClust = filoBranch(frame).reconstructInfo.CandMaskPreCluster;
idx = find(preClust);
[preClustY,preClustX] = ind2sub(size(preClust),idx);
c = brewermap(2,'dark2');
scatter(preClustX,preClustY,5,c(2,:),'filled');

if ip.Results.plotScaleBar
    plotScaleBar(pixels,pixels/10,'Color',textColor);
end

countFig = countFig+1;
%% Show Clustering
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
spy(seedMask,'b');
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
if ip.Results.writeTitles
    text(nx/10, 10,'Cluster Linear Candidates', fontText{:});
end
clusterLinks = filoBranch(frame).reconstructInfo.clusterlinks;

preClust = filoBranch(frame).reconstructInfo.CandMaskPreCluster;
idx= find(preClust);
[preClustY,preClustX] = ind2sub(size(preClust),idx);
scatter(preClustX,preClustY,5,c(2,:),'filled');
spy(clusterLinks,'k');

if ip.Results.plotScaleBar
    plotScaleBar(pixels,pixels/10,'Color',textColor);
end

countFig = countFig+1;
%% Candidates Post Clustering
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
if ip.Results.writeTitles
    text(nx/10, 10,' Linear Candidates Clustered', fontText{:});
end
spy(seedMask,'b');
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
postClust= filoBranch(frame).reconstructInfo.CandMaskPostCluster;
idx = find(postClust);
[postClustY,postClustX] = ind2sub(size(postClust),idx);
scatter(postClustX,postClustY,5,c(2,:),'filled');
pixels = round(10/pixSizeMic);

if ip.Results.plotScaleBar
    plotScaleBar(pixels,pixels/10,'Color',textColor);
end

countFig = countFig+1;
%% Iterate over reconstruction
if isfield(filoBranch(frame).reconstructInfo,'output');
    for iReconst = 1:numel(filoBranch(frame).reconstructInfo.output)
        
        hSet(countFig).h = setFigure(nx,ny) ;
        imshow(-img,[])
        hold on
        if ip.Results.writeTitles
            text(nx/10, 10,'Link Candidates', fontText{:});
        end
        
        idx = find(postClust);
        [postClustY,postClustX] = ind2sub(size(postClust),idx);
        scatter(postClustX,postClustY,5,c(2,:),'filled');
        spy(seedMask,'b');
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        
        links= filoBranch(frame).reconstructInfo.output{iReconst}.links;
        spy(links,'k',5);
      
        if ip.Results.plotScaleBar
            plotScaleBar(pixels,pixels/10,'Color',textColor);
        end
       
        countFig  = countFig+1;
        %% Add to VeilStem
        hSet(countFig).h = setFigure(nx,ny) ;
        
        imshow(-img,[])
        hold on
        
        spy(seedMask,'b');
     
        if ip.Results.writeTitles
            text(nx/10, 10,'Add to VeilStem', fontText{:});
        end
        bodyAdd =  filoBranch(frame).reconstructInfo.output{iReconst}.candFiloAdded.Body ;
        
        [yB,xB] = ind2sub(size(bodyAdd),find(bodyAdd));
        
        aqua = [0.0039  ,  0.4264 ,   0.3848];
        scatter(xB,yB,8,c(2,:),'filled');
        
        if ip.Results.plotScaleBar
            plotScaleBar(pixels,pixels/10,'Color',textColor);
        end
        
        countFig = countFig+1;
        %% 'Add Branch'
        hSet(countFig).h = setFigure(nx,ny);
        imshow(-img,[])
        hold on
        
        spy(seedMask,'b');
        scatter(xB,yB,5,'b','filled');
        
        if ip.Results.writeTitles
            text(nx/10, 10,'Add Branch', fontText{:});
        end
        branchAdd =  filoBranch(frame).reconstructInfo.output{iReconst}.candFiloAdded.Branch ;
        if sum(branchAdd(:))>0
            [yBr,xBr] = ind2sub(size(branchAdd),find(branchAdd));
           
            scatter(xBr,yBr,8,c(2,:),'filled');
        end
        
        if ip.Results.plotScaleBar
            plotScaleBar(pixels,pixels/10,'Color',textColor);
        end
        
        countFig  = countFig +1;
        %% Add End-on Attachment
        hSet(countFig).h = setFigure(nx,ny);
        imshow(-img,[])
        hold on
        
        spy(seedMask,'b');
        scatter(xB,yB,5,'b','filled');
        
        if ip.Results.writeTitles
            text(nx/10, 10,'Add End-On Attachment', fontText{:});
        end
        endOn =  filoBranch(frame).reconstructInfo.output{iReconst}.candFiloAdded.EndOn ;
        
        if sum(endOn(:))>0
            [yE,xE] = ind2sub(size(endOn),find(endOn));
            scatter(xE,yE,8,c(2,:),'filled');
        end
        
        if ip.Results.plotScaleBar
            plotScaleBar(pixels,pixels/10,'Color',textColor);
        end
        
        countFig = countFig+1;        
        %%
        if iReconst == numel(filoBranch(frame).reconstructInfo.output);
            
            title = 'End Reconstruction';
            
        else
            title = 'New Seed';
        end
        hSet(countFig).h = setFigure(nx,ny) ;
        imshow(-img,[])
        hold on
        
        if ip.Results.writeTitles
            text(nx/10, 10,title, fontText{:});
        end
        
        seedMask =  filoBranch(frame).reconstructInfo.seedMask{iReconst+1} ;
        spy(seedMask,'b')
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        
        if ip.Results.plotScaleBar
            plotScaleBar(pixels,pixels/10,'Color',textColor);
        end
   
        countFig = countFig +1;
    end % iReconst
end % isfield

filoInfo = filoBranch(frame).filoInfo;


if isfield(filoInfo,'Ext_params'); 
    
%     %%   Filopodia By Branch Group
%     
%     hSet(countFig).h  = setFigure(nx,ny);
%     imshow(-img,[]);
%     hold on
%     
%     GCAVisualsPlotFilopodiaPerBranchGroup(filoInfo,[ny,nx]);
%     cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);
%     
%     if ip.Results.writeTitles
%         text(nx/10,10,'Color-Coded By Branch Group', fontText{:});
%     end
%     
%     countFig  = countFig+1;
%%    
%     %% Plot Individual: For Validation : Random Permutation of Colors
%     hSet(countFig).h = setFigure(nx,ny);
%     imshow(-img,[]);
%     hold on
%     cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);
%     
%     [filoFilterSet,filoParams] = GCACreateFilopodiaFilterSetWithEmbedResidTest(filoBranch,'Validation_NoEmbed');
%     
%     filoFilterSetC = filoFilterSet(frame);
%     
    % Filter the filopodia
%     filoInfoG = filoInfo(filoFilterSetC{1}(:,1));
%     n = length(filoInfoG);
%     
%     %c = linspecer(n);
%     c = brewermap(n,'spectral');
%     
%     idxRand = randperm(n);
%     c = c(idxRand,:);
%     
%     % plot for each filo
%     for ifilo = 1:length(filoInfoG)
%         filoInfoC = filoInfoG(ifilo);
%         GCAVisualsMakeOverlaysFilopodia(filoInfoC,[ny,nx],0,1,c(ifilo,:),0);
%         clear filoInfoC
%     end
%     
%     % plot filo documented but did not pass the criteria in black
%     filoInfoB = filoInfo(~filoFilterSetC{1}(:,1));
%     % GCAVisualsMakeOverlaysFilopodia(filoInfoB,[ny,nx],0,1,'k',0,1);
%     
%     
%     if ip.Results.writeTitles
%         text(nx/10,10,'Color-Coded Individual Segment', fontText{:});
%     end
%     countFig = countFig+1;
%     
    %% ColorCode by Length Using either the brewermap or the jet map -
 
    %cMap{1} = flip(brewermap(128,'spectral'),1);
    
    cMap{1} = jet(128);
    
    filoBranchC = filoBranch(frame);
    
    [filoFilterSet,filoParams] = GCACreateFilopodiaFilterSetWithEmbedResidTest(filoBranch,'Validation_NoEmbed');
    
    filoFilterSetC = filoFilterSet(frame);
    
    filoLengths = GCAAnalysisExtract_filoLength(filoBranchC,filoFilterSetC,'umPerPixel',pixSizeMic);
  
    plotValues = filoLengths{1};
    cMapLimits(1) = 0;
    cMapLimits(2) = 10;
  
    for iMap = 1:1
        
        hSet(countFig).h = setFigure(nx,ny);
        imshow(-img,[]);
        hold on
        cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);
        
        
        imgSize = size(img);
        GCAVisualsFilopodiaMeasurementOverlays(filoInfo,imgSize,...
            'filoFilterSet',filoFilterSetC{1},'plotValues',plotValues,...
            'justExt',1,'cMapLimits',cMapLimits,'ColorByValue',true, ...
            'plotText',false,'colorMap',cMap{iMap},'ExtraColor',[]);
       
         if ip.Results.plotScaleBar
            plotScaleBar(pixels,pixels/10,'Color',textColor);
        end
        
        countFig = countFig+1;
        
    end % iMap
    
else
    filoFilterSet = [];
    filoParams = [];
end

end
