function [ hSet ] = GCATroubleShootMakeMovieOfReconstruct( filoBranch,veilStem,frame,pixSizeMic,imDir)
% Small function that makes a cool movie of the different steps in the reconstruction
% likely also will be very helpful for troubleshooting

% if otherImgDir == 0;
% filename = analInfo(frame).imgPointer;
% textColor = [ 1 1 1];
% else
listOfImages = searchFiles('.tif',[],imDir);
filename = [char(listOfImages(frame,2)) filesep char(listOfImages(frame,1))];
textColor = [ 0 0 0 ];
% end
fontText =  {'FontName','Arial','FontSize',14,'FontName','Arial','color',textColor};
%%%% START %%%%
img = double(imread(filename));
% img = [img img];
[ny,nx] = size(img);

%% Initiate 
countFig = 1; 
%% 00 -original image

hSet(countFig).h =setFigure(nx,ny);
imshow(-img,[]) ;
hold on
text(nx/10, 10,'Original Image ', fontText{:});
pixels = 10/pixSizeMic;
plotScaleBar(pixels,pixels/10,'Color',textColor);
countFig = countFig+1; 
%% 01 - LOCAL THRESH %%%%
% h = setFigure(nx,ny);
% imshow(-img,[])
% hold on
%
% localThresh = analInfo(frame).bodyEst.beforeErod;
%
% cellfun(@(x) plot(x(:,2),x(:,1),'c'),localThresh);
%
% text(nx/10, 10,'Step 01: Local Thresholding', 'color',textColor);
% pixels = round(10/pixSizeMic);
%      plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
%
%    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
%             [saveDir filesep 'ReconstructMovie01.png']);
%         saveas(h,[saveDir filesep 'ReconstructMovie01.eps'],'psc2');

%         close gcf

%% 02 - EROD FILO %%%%
% h = setFigure(nx,ny) ;
% imshow(-img,[])
% hold on
% erosion =veilStem(frame).bodyEst.erodForBody;
% cellfun(@(x) plot(x(:,2),x(:,1),'b'),erosion)
% text(nx/10, 20,{'Estimate Larger-Scale'; 'Veil/Stem Pieces'}, fontText{:});
% pixels = round(10/pixSizeMic);
% plotScaleBar(pixels,pixels/20,'Color',textColor);
% 
% %  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
% %          [saveDir filesep 'ReconstructMovie02.png']);
% saveas(h,[saveDir filesep 'ReconstructMovie02.png']);
% saveas(h,[saveDir filesep 'ReconstructMovie02.eps'],'psc2');
% saveas(h,[saveDir filesep 'ReconstructMovie02.fig']);
% 
% close gcf

%% 03 BACKBONE SAVE %%%%
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
backbone =  veilStem(frame).backbone;
% extra = zeros(size(img));
% backbone = [backbone extra];
spy(backbone,'r')
%cellfun(@(x) plot(x(:,2),x(:,1),'--b'),erosion)
text(nx/10, 20,{'Estimate Backbone:' ; 'Large Scale Ridges'}, fontText{:});
pixels = round(10/pixSizeMic);
plotScaleBar(pixels,pixels/20,'Color',textColor);

countFig = countFig+1; 
%% 04 Final Body Mask %%%%
hSet(countFig).h = setFigure(nx,ny);
imshow(-img,[]) ;
hold on
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
text(nx/10, 10,'Veil/Stem Estimation Complete', fontText{:});
pixels = round(10/pixSizeMic);
plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
countFig = countFig+1; 
%% 05 Overlay Ridges %%%%
% h = setFigure(nx,ny);
% imshow(-img,[]) ;
% hold on
% candRidges = filoBranch(frame).filterInfo.ThreshNMS;
% spy(candRidges ,'m');
% cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
% text(nx/10, 10,'Detect Small Scale Ridges', fontText{:});
% pixels = round(10/pixSizeMic);
% %  plotScaleBar(pixels,pixels/20,'Color',textColor);
% % print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
% %         [saveDir filesep 'ReconstructMovie05.png']);
% saveas(h,[saveDir filesep 'ReconstructMovie05.png'] );
% saveas(h,[saveDir filesep 'ReconstructMovie05.eps'],'psc2');

% close gcf

%% 06 Show Seed
hSet(countFig).h = setFigure(nx,ny);
imshow(-img,[]) ;
hold on


text(nx/10, 10,'Get Seed For Reconstruction', fontText{:});
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
seedMask = filoBranch(frame).reconstructInfo.seedMask{1};
spy(seedMask,'b',5);
pixels = round(10/pixSizeMic);
 plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
countFig = countFig+1; 
%% 07 Show Candidates
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
text(nx/10, 10,'Get Candidates', fontText{:});
spy(seedMask,'b');
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
preClust = filoBranch(frame).reconstructInfo.CandMaskPreCluster;

spy(preClust,'m',5);
pixels = round(10/pixSizeMic);
plotScaleBar(pixels,pixels/20,'Color',textColor);
countFig = countFig+1; 

%% 08 Show Clustering
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
spy(seedMask,'b');
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
text(nx/10, 10,'Cluster Linear Candidates', fontText{:});
clusterLinks = filoBranch(frame).reconstructInfo.clusterlinks;

preClust = filoBranch(frame).reconstructInfo.CandMaskPreCluster;
spy(preClust,'m');
spy(clusterLinks,'y');
pixels = round(10/pixSizeMic);
plotScaleBar(pixels,pixels/20,'Color',textColor);

countFig = countFig+1; 
%% 09 Candidates Post Clustering
hSet(countFig).h = setFigure(nx,ny) ;
imshow(-img,[])
hold on
text(nx/10, 10,' Linear Candidates Clustered', fontText{:});
spy(seedMask,'b');
bodyFinal = veilStem(frame).finalMask;
edgeYX = bwboundaries(bodyFinal);
cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
postClust= filoBranch(frame).reconstructInfo.CandMaskPostCluster;
spy(postClust,'m');
pixels = round(10/pixSizeMic);
plotScaleBar(pixels,pixels/10,'Color',textColor);
countFig = countFig+1; 

%% Iterate over reconstruction
if isfield(filoBranch(frame).reconstructInfo,'output');
    for iReconst = 1:numel(filoBranch(frame).reconstructInfo.output)
        
        hSet(countFig).h = setFigure(nx,ny) ;
        imshow(-img,[])
        hold on
        text(nx/10, 10,'Link Candidates', fontText{:});
        spy(postClust,'m');
        spy(seedMask,'b');
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        
        links= filoBranch(frame).reconstructInfo.output{iReconst}.links;
        spy(links,'y',5);
        pixels = round(10/pixSizeMic);
        plotScaleBar(pixels,pixels/10,'Color',textColor); 
        countFig  = countFig+1; 
        %%
        hSet(countFig).h = setFigure(nx,ny) ;
        
        imshow(-img,[])
        hold on
        spy(postClust,'m');
        spy(seedMask,'b');
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        
        spy(links,'b');
        text(nx/10, 10,'Add to Body', fontText{:});
        bodyAdd =  filoBranch(frame).reconstructInfo.output{iReconst}.candFiloAdded.Body ;
        spy(bodyAdd,'g',5)
        pixels = round(10/pixSizeMic);
         plotScaleBar(pixels,pixels/10,'Color',textColor);
        countFig = countFig+1; 
        %%
        hSet(countFig).h = setFigure(nx,ny); 
        imshow(-img,[])
        hold on
        spy(postClust,'m');
        spy(seedMask,'b');
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        
        spy(links,'b');
        spy(bodyAdd,'b')
        text(nx/10, 10,'Add Branch', fontText{:});
        branchAdd =  filoBranch(frame).reconstructInfo.output{iReconst}.candFiloAdded.Branch ;
        spy(branchAdd,'g',5)
        pixels = round(10/pixSizeMic);
          plotScaleBar(pixels,pixels/10,'Color',textColor);
        
        countFig  = countFig +1; 
        %%
        hSet(countFig).h = setFigure(nx,ny); 
        imshow(-img,[])
        hold on
        
        spy(postClust,'b');
        spy(seedMask,'b');
        
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        spy(links,'b');
        spy(bodyAdd,'b')
        spy(branchAdd,'b');
        text(nx/10, 10,'Add End-On Attachment', fontText{:});
        endOn =  filoBranch(frame).reconstructInfo.output{iReconst}.candFiloAdded.EndOn ;
        spy(endOn,'g',5)
        pixels = round(10/pixSizeMic);
        plotScaleBar(pixels,pixels/10,'Color',textColor);
       
        countFig = countFig+1; 
        
        %%
        if iReconst == numel(filoBranch(frame).reconstructInfo.output);
            
            title = 'End Reconstruction';
            
        else
            title = 'New Seed';
        end
        hSet(countFig) = setFigure(nx,ny) ;
        imshow(-img,[])
        hold on
        text(nx/10, 10,title, fontText{:});
        seedMask =  filoBranch(frame).reconstructInfo.seedMask{iReconst+1} ;
        spy(seedMask,'b')
        bodyFinal = veilStem(frame).finalMask;
        edgeYX = bwboundaries(bodyFinal);
        cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX);
        
        pixels = round(10/pixSizeMic);
        plotScaleBar(pixels,pixels/10,textColor);
        
       
        countFig = countFig+1; 
    end % iReconst
end
%% 
% %% plot fits 1st, 2nd, and higher order branches
%
%   h  = setFigure(nx,ny);
% %figure;
%
%  imshow(-img,[]);
%  hold on
%
%  cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX);
%  filoInfoSingle = filoInfo(type==0);
%  plotfilosIntAndExt(filoInfoSingle,imgSize,1,1,'g',0);
%
%  text(nx/10,10,'High Confidence Single Filopodia Attached to Neurite Body', 'color',textColor);
%  pixels = round(10/pixSizeMic);
%    plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
%
%
%
%
%
%         test = vertcat(filoInfoSingle(:).Ext_exitFlag);
%         %c(iType) = c;
%         idx = test>=1;
%
%         filoInfoFilt = filoInfo(idx);
%
%         % filter out any that might have passed the exitflag criteria but NOT
%         % gave a number for the fit ==0  % maybe flag above later...
%         test2 = vertcat(filoInfoFilt(:).Ext_endpointCoordFitPix);
%         idx2 = ~isnan(test2);
%         filoInfoFilt = filoInfoFilt(idx2);
%    value =  nanmean(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic;
%
%
%   text(nx/10,30,['Mean Length ' num2str(value,2) ' um'], 'color',textColor);
%    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
%             [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png'])
%  clear filoInfoFilt
%         imageNum = imageNum+1;
%         close gcf
%%
%  h  = setFigure(nx,ny);
% %figure;
%
%  imshow(-img,[]);
%  hold on
%
%  cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX);
%  filoInfoBranchStem = filoInfo(type==1);
%  plotfilosIntAndExt(filoInfoBranchStem,imgSize,1,1,'y',1);
%
%  text(nx/10,10,'High Confidence Branch Stem Attached to Neurite Body', 'color',textColor);
%  pixels = round(10/pixSizeMic);
%   plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
%
%    test = vertcat(filoInfoBranchStem(:).Ext_exitFlag);
%         %c(iType) = c;
%         idx = test>=1;
%
%         filoInfoFilt = filoInfoBranchStem(idx);
%
%         % filter out any that might have passed the exitflag criteria but NOT
%         % gave a number for the fit ==0  % maybe flag above later...
%         test2 = vertcat(filoInfoFilt(:).Ext_endpointCoordFitPix);
%         idx2 = ~isnan(test2);
%         filoInfoFilt = filoInfoFilt(idx2);
%    value =  nanmean(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic;
%     valueMax = nanmax(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic;
%   text(nx/10,30,['Mean Length ' num2str(value,2) ' um'], 'color',textColor);
%    text(nx/10,50,['Max Length ' num2str(valueMax,2) ' um'],'color',textColor);
%    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
%             [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png'])
%
%         close gcf
%         imageNum = imageNum+1;
%         clear filoInfoFilt
%%
% h  = setFigure(nx,ny);
% %figure;
%
%  imshow(-img,[]);
%  hold on
%
%  cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX);
%  filoInfoHigherOrder = filoInfo(type==2);
%  plotfilosIntAndExt(filoInfoHigherOrder,imgSize,1,1,'m',1);
%
%  text(nx/10,10,'High Confidence 1st Order Branches', 'color',textColor);
%  pixels = round(10/pixSizeMic);
%    plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
%  test = vertcat(filoInfoHigherOrder(:).Ext_exitFlag);
%         %c(iType) = c;
%         idx = test>=1;
%
%         filoInfoFilt = filoInfoHigherOrder(idx);
%
%         % filter out any that might have passed the exitflag criteria but NOT
%         % gave a number for the fit ==0  % maybe flag above later...
%         test2 = vertcat(filoInfoFilt(:).Ext_endpointCoordFitPix);
%         idx2 = ~isnan(test2);
%         filoInfoFilt = filoInfoFilt(idx2);
%    value =  nanmean(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic;
%     valueMax = nanmax(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic;
%   text(nx/10,30,['Mean Length ' num2str(value,2) ' um'], 'color',textColor);
%   text(nx/10,50,['Max Length ' num2str(valueMax,2) ' um'],'color',textColor);
%    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
%             [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png'])
%         close gcf

filoInfo = filoBranch(frame).filoInfo;

%%   Filopodia By Branch Group

hSet(countFig)  = setFigure(nx,ny);
imshow(-img,[]);
hold on

GCAVisualsPlotFilopodiaPerBranchGroup(filoInfo,[ny,nx]);
cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);

text(nx/10,10,'Color-Coded By Branch Group', fontText{:});
countFig  = countFig+1; 

%% Plot Individual
hSet(countFig).h = setFigure(nx,ny);
imshow(-img,[]);
hold on
cellfun(@(x) plot(x(:,2),x(:,1),'k'),edgeYX);

n = length(filoInfo);
c = linspecer(n);
idxRand = randperm(n);
c = c(idxRand,:);

for ifilo = 1:length(filoInfo)
    filoInfoC = filoInfo(ifilo);
    GCAVisualsMakeOverlaysFilopodia(filoInfoC,[ny,nx],1,1,c(ifilo,:),0);
    clear filoInfoC
end


text(nx/10,10,'Color-Coded Individual Segment', fontText{:});

countFig  = countFig+1; 

% cd(saveDir)
% execute = 'mencoder mf://*.png -mf w=800:h=600:fps=0.5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
% system(execute);





end

