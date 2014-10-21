function [ output_args ] = makeMovieOfReconstruct( analInfo,frame,pixSizeMic,saveDir,imDir)
% Small function that makes a cool movie of the different steps in the reconstruction 
% likely also will be very helpful for troubleshooting

% if otherImgDir == 0; 
% filename = analInfo(frame).imgPointer; 
% textColor = [ 1 1 1]; 
% else 
    listOfImages = searchFiles('.tif',[],imDir); 
    filename = [char(listOfImages(frame,2)) filesep char(listOfImages(frame,1))];
    textColor = [ 0  0 0 ]; 
% end 
 
%%%% START %%%%     
img = double(imread(filename)); 
% img = [img img]; 
[ny nx] = size(img); 
imgSize = [ny nx]; 
zoom = 1; 

%% 00 -original image
h=setFigure(nx,ny); 
imshow(-img,[]) ; 
hold on 
text(nx/10, 10,'Original Image', 'color',textColor);
pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie00.png']); 
       
saveas(h,[saveDir filesep 'ReconstructMovie00.eps'],'psc2')
 close gcf
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
h = setFigure(nx,ny) ;
imshow(-img,[]) 
hold on 
  erosion =analInfo(frame).bodyEst.erodForBody; 
cellfun(@(x) plot(x(:,2),x(:,1),'b'),erosion)
text(nx/10, 10,'Step 02: Erosion', 'color',textColor);
pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie02.png']); 
          saveas(h,[saveDir filesep 'ReconstructMovie02.eps'],'psc2'); 
        
        close gcf
      
%% 03 BACKBONE SAVE %%%%
h = setFigure(nx,ny) ;
imshow(-img,[]) 
hold on 
backbone =  analInfo(frame).bodyEst.backbone; 
% extra = zeros(size(img)); 
% backbone = [backbone extra]; 
spy(backbone,'r')
cellfun(@(x) plot(x(:,2),x(:,1),'--b'),erosion)
text(nx/10, 10,'Step 03: Backbone Estimate', 'color',textColor);
pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie03.png']); 
        saveas(h,[saveDir filesep 'ReconstructMovie03.eps'],'psc2'); 
close gcf 
%% 04 Final Body Mask %%%% 
  h = setFigure(nx,ny); 
  imshow(-img,[]) ; 
  hold on 
  bodyFinal = analInfo(frame).masks.neuriteEdge; 
  edgeYX = bwboundaries(bodyFinal); 
  cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX); 
  text(nx/10, 10,'Step 04: Final Mask', 'color',textColor);
    pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie04.png']);
        saveas(h,[saveDir filesep 'ReconstructMovie04.eps'],'psc2'); 
 close gcf 
%% 05 Overlay Ridges %%%% 
 h = setFigure(nx,ny); 
  imshow(-img,[]) ; 
  hold on 
  candRidges = analInfo(frame).filterInfo.ThreshNMS;
  spy(candRidges ,'m');
  cellfun(@(x) plot(x(:,2),x(:,1),'b'),edgeYX); 
   text(nx/10, 10,'Step 05: Ridges!', 'color',textColor);
    pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie05.png']);
          saveas(h,[saveDir filesep 'ReconstructMovie05.eps'],'psc2'); 
        
 close gcf
  
  %% 06 Show Seed 
  h = setFigure(nx,ny); 
  imshow(-img,[]) ; 
  hold on 
 
  
  text(nx/10, 10,'Step 06: Get Seed', 'color',textColor);
  seedMask = analInfo(frame).reconstructInfo.seedMask{1}; 
  spy(seedMask,'b',5); 
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie06.png']);
            saveas(h,[saveDir filesep 'ReconstructMovie06.eps'],'psc2'); 
  close gcf 
  %% 07 Show Candidates
  h = setFigure(nx,ny) ;
  imshow(-img,[]) 
  hold on 
  text(nx/10, 10,'Step 07: Get Candidates', 'color',textColor);
  spy(seedMask,'b'); 
  preClust = analInfo(frame).reconstructInfo.CandMaskPreCluster; 
  
  spy(preClust,'m',5); 
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie07.png']);
            saveas(h,[saveDir filesep 'ReconstructMovie07.eps'],'psc2'); 
        saveas(h,[saveDir filesep 'ReconstructMovie07.eps'],'psc2'); 
        
%% 08 Show Clustering 
h = setFigure(nx,ny) ;
imshow(-img,[]) 
  hold on 
    spy(seedMask,'b'); 
  text(nx/10, 10,'Step 08: Cluster Linear Candidates', 'color',textColor);
  clusterLinks = analInfo(frame).reconstructInfo.clusterlinks; 
  
  preClust = analInfo(frame).reconstructInfo.CandMaskPreCluster; 
  spy(preClust,'m'); 
  spy(clusterLinks,'y'); 
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie08.png']);
        saveas(h,[saveDir filesep 'ReconstructMovie08.eps'],'psc2'); 
close gcf
%% 09 Candidates Post Clustering 
h = setFigure(nx,ny) ;
imshow(-img,[]) 
  hold on 
  text(nx/10, 10,' Linear Candidates Clustered', 'color',textColor);
   spy(seedMask,'b'); 
  postClust= analInfo(frame).reconstructInfo.CandMaskPostCluster; 
  spy(postClust,'m'); 
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie09.png']);
          saveas(h,[saveDir filesep 'ReconstructMovie09.eps'],'psc2'); 
close gcf
%% Iterate over reconstruction 
imageNum = 10; 
for iReconst = 1:numel(analInfo(frame).reconstructInfo.output)

h = setFigure(nx,ny) ;
imshow(-img,[]) 
  hold on 
  text(nx/10, 10,'Step 09: Link Candidates', 'color',textColor);
 spy(postClust,'m'); 
  spy(seedMask,'b'); 
 
 links= analInfo(frame).reconstructInfo.output{iReconst}.links; 
  spy(links,'y',5); 
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png']);
          saveas(h,[saveDir filesep 'ReconstructMovie' num2str(imageNum) '.eps'],'psc2'); 
        imageNum = imageNum +1; 
close gcf

%% 
h = setFigure(nx,ny) ;

imshow(-img,[]) 
  hold on 
   spy(postClust,'m'); 
  spy(seedMask,'b'); 
 
  spy(links,'b'); 
  text(nx/10, 10,'Add to Body', 'color',textColor);
  bodyAdd =  analInfo(frame).reconstructInfo.output{iReconst}.candFiloAdded.Body ; 
  spy(bodyAdd,'g',5)
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
 print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png']);
        imageNum = imageNum+1; 
%%
imshow(-img,[]) 
  hold on 
   spy(postClust,'m');
  spy(seedMask,'b'); 
 
  spy(links,'b'); 
  spy(bodyAdd,'b')
  text(nx/10, 10,'Add Branch', 'color',textColor);
 branchAdd =  analInfo(frame).reconstructInfo.output{iReconst}.candFiloAdded.Branch ; 
  spy(branchAdd,'g',5)
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
 print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png']);  
        imageNum = imageNum+1; 
%%        
  imshow(-img,[]) 
  hold on 
  
  spy(postClust,'b'); 
  spy(seedMask,'b'); 
  spy(links,'b');
   spy(bodyAdd,'b')
   spy(branchAdd,'b'); 
  text(nx/10, 10,'End On Attachment', 'color',textColor);
  endOn =  analInfo(frame).reconstructInfo.output{iReconst}.candFiloAdded.EndOn ; 
  spy(endOn,'g',5)
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
 print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) ' .png']);  
          saveas(h,[saveDir filesep 'ReconstructMovie' num2str(imageNum) ' .eps'],'psc2');  
  imageNum = imageNum+1; 
  close gcf
%%
if iReconst == numel(analInfo(frame).reconstructInfo.output); 

    title = 'End Reconstruction';
    
else 
    title = 'New Seed';
end 
h = setFigure(nx,ny) ;
imshow(-img,[]) 
  hold on 
  text(nx/10, 10,title, 'color',textColor);
  seedMask =  analInfo(frame).reconstructInfo.seedMask{iReconst+1} ; 
  spy(seedMask,'b')
 
  pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png']);
        saveas(h,[saveDir filesep 'ReconstructMovie' num2str(imageNum) ' .eps'],'psc2');
        imageNum = imageNum+1; 
       
       close gcf
end % iReconst
%% sanity check show crosses 


h  = setFigure(nx,ny); 
 filoInfo = analInfo(frame).filoInfo; 
  filoInfo = replaceBadFitsWithOriginalLength(filoInfo,imgSize);
filoInfo = addFiloFitCoords(filoInfo); 
%c = ['r' 'b' 'g' 'y' 'c' 'm' 'k','r']; 
%c = c'; 
 c = colormap(lines(numel(filoInfo)));
%  sums = arrayfun(@(x) sum(c(x,:)),1:length(c(:,1))); 
%  idxBad = find(sums(:)==0.75); % take out black
%  x = repmat([0.5 0 0.5],length(idxBad),1); 
% c(idxBad,:) = x; 
 imshow(-img,[]); 
 hold on 
   
    for i = 1:numel(filoInfo)
%      pixIndices = filoInfo(i).Ext_pixIndices;
%             idxEnd = find(pixIndices == filoInfo(i).Ext_pixIndicesendpointCoordFitPix;
%             pixIndicesPlot = pixIndices(1:idxEnd);   
%         ind2sub(imgSize)
    plot(filoInfo(i).Ext_toPlotXY(:,1),filoInfo(i).Ext_toPlotXY(:,2),'color',c(i,:))
  %  scatter(filoInfo(i).Ext_endpointCoordFitXY(:,1),filoInfo(i).Ext_endpointCoordFitXY(:,2),50,c(i,:),'filled'); 
    end 

pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
text(nx/10,10,'Color-Coded by Filopodia Number','color',textColor); 
 print(h, '-depsc', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.eps']);
imageNum = imageNum+1; 

close gcf
%% show BranchGroups

 

%  connIdx = arrayfun(@(x) ~isempty(x.conIdx),filoInfo);
%  branches = filoInfo(connIdx); 
 type = arrayfun(@(x) x.type,filoInfo); 
% % filoInfoBranch = filoInfo(type==1); 
%  
% %  sums = arrayfun(@(x) sum(c(x,:)),1:length(c(:,1))); 
% %  idxBlack = find(sums(:)==.75); 
% %   x = repmat([0.5 0 0.5],length(idxBlack),1); 
% % c(idxBlack,:) = x; 
% 
% filoInfoBranch = filoInfo(type ~=0); 
% groupIDAll =vertcat(filoInfoBranch(:).groupCount); 
% groupID = unique(groupIDAll); 
%  c = colormap(lines(length(groupID)));
% c = ['r' 'b' 'g' 'y' 'c' 'm' 'k','r']; 
% c = c'; 
% h  = setFigure(nx,ny); 
% 
%  imshow(-img,[]); 
%  hold on 
%  
%  text(nx/10,10,'Color-Coded by BranchGroup','color',textColor); 
%  
%  for iGroup = 1:length(groupID)
%    %grpIdx = find(groupIDAll  == groupID(iGroup)); 
%    filoInfoGroupC = filoInfoBranch(groupIDAll==groupID(iGroup)); 
%    
%   
%    for i = 1:numel(filoInfoGroupC)
%         plot(filoInfoGroupC(i).Ext_toPlotXY(:,1),filoInfoGroupC(i).Ext_toPlotXY(:,2),'color',c(iGroup,:))
%        scatter(filoInfoGroupC(i).Ext_endpointCoordFitXY,filoInfoGroupC(i).Ext_endpointCoordFitXY,20,c(iGroup,:));
%    end 
%    
%    
%  end 
%    
%     
%  
 
 
%  for i = 1:numel(filoInfoBranch)
%  plotfilosIntAndExt(filoInfoBranch(i),[ny,nx],0,1,c(i,:)); 
%  conIdx = filoInfoBranch(i).conIdx; 
%  plotfilosIntAndExt(filoInfo(conIdx),[ny,nx],0,1,c(i,:));
%  end 
%  % this is a stupid way to save the branch data need to think of
%  something more clever. 
% get max type 
% find max type 
 
% 
%     filoInfoBranchMain = filoInfo(type==1); % first just get the main branches 
%     
%    
%  for i = 1:numel(filoInfoBranchBranchMain);  
%      % get xyCoords 
%      plot(filoInfoBranchMain(i).Ext_toPlotXY(:,1),filoInfoBranchMain(i).Ext_toPlotXY(:,2),'color',c(i,:)); 
%      
%      % get connectivity of that main branch % this is the ID of all
%      % connecting filopodia processes
%      conIdx = filoInfoBranch(i).conIdx;
%      % plot all branches in same color
%      for j = 1:numel(conIdx)
%          plot(filoInfo(conIdx(j)).Ext_toPlotXY(:,1),filoInfo(conIdx(j)).Ext_toPlotXY(:,2),'color',c(i,:)); 
%      end 
%      % scatter endpoints 
%       xycoordsEndFitMainBranch  = filoInfoBranch(i).Ext_endpointCoordFitXY;
%       scatter(xycoordsEndFitMainBranch(:,1),xycoordsEndFitMainBranch(:,2),20,c(i,:),'filled'); 
%       xycoordsEndFitExtraBranch = vertcat(filoInfo(conIdx).Ext_endpointCoordFitXY); 
%       scatter(xycoordsEndFitExtraBranch(:,1),xycoordsEndFitExtraBranch(:,2),20,c(i,:),'filled'); 
%       xyCoordsCon = vertcat(filoInfo(conIdx).conXYCoords); 
%       if ~isempty(xyCoordsCon)
%       scatter(xyCoordsCon(:,1),xyCoordsCon(:,2),20,c(i,:),'filled'); 
%       end 
%     
%     clear  conIdx
%     
%  end 
%  clear filoInfoBranch
%  iBranchType = iBranchType +1 ; 
% end 
 
 
%  
%  hold on 
%      branches = vertcat(filoInfo(:).conXYCoords);
%      idxBranches = sub2ind(imgSize,branches(:,2),branches(:,1)); 
%      numBranchPoints = length(unique(idxBranches)); 
% if ~isempty(branches)
%   scatter(branches(:,1),branches(:,2),20,'y','filled'); 
% end 
% pixels = round(10/pixSizeMic); 
%      plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
%      
% text(nx/10,20,[num2str(numBranchPoints) ' Branch Junctions'],'color',textColor); 
%  print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
%             [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png'])
%         imageNum = imageNum+1;
%        close gcf 
 %% plot fits 
 h  = setFigure(nx,ny); 

 imshow(-img,[]); 
 hold on
 text(nx/10,10,'Higher Confidence Length Measurements', 'color',textColor);
 pixels = round(10/pixSizeMic); 
     plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
     
  cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX); 
  hold on 
    plotfilosIntAndExt(filoInfo,[ny,nx],1,1,[]);
    print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png'])
            saveas(h,[saveDir filesep 'ReconstructMovie' num2str(imageNum) '.eps'],'psc2'); 
        imageNum = imageNum+1; 
        close gcf
%% plot fits 1st, 2nd, and higher order branches         
 
  h  = setFigure(nx,ny); 
%figure; 

 imshow(-img,[]);   
 hold on 
 
 cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX); 
 filoInfoSingle = filoInfo(type==0); 
 plotfilosIntAndExt(filoInfoSingle,imgSize,1,1,'g',0); 
 
 text(nx/10,10,'High Confidence Single Filopodia Attached to Neurite Body', 'color',textColor);
 pixels = round(10/pixSizeMic); 
   plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
 
 
 
    
   
        test = vertcat(filoInfoSingle(:).Ext_exitFlag);
        %c(iType) = c;
        idx = test>=1;
        
        filoInfoFilt = filoInfo(idx);
        
        % filter out any that might have passed the exitflag criteria but NOT
        % gave a number for the fit ==0  % maybe flag above later...
        test2 = vertcat(filoInfoFilt(:).Ext_endpointCoordFitPix);
        idx2 = ~isnan(test2);
        filoInfoFilt = filoInfoFilt(idx2);
   value =  nanmean(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic; 
    
  
  text(nx/10,30,['Mean Length ' num2str(value,2) ' um'], 'color',textColor);     
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png'])  
 clear filoInfoFilt
        imageNum = imageNum+1; 
        close gcf
 %%       
 h  = setFigure(nx,ny); 
%figure; 

 imshow(-img,[]);   
 hold on 
 
 cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX); 
 filoInfoBranchStem = filoInfo(type==1); 
 plotfilosIntAndExt(filoInfoBranchStem,imgSize,1,1,'y',1); 
 
 text(nx/10,10,'High Confidence Branch Stem Attached to Neurite Body', 'color',textColor);
 pixels = round(10/pixSizeMic);         
  plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
 
   test = vertcat(filoInfoBranchStem(:).Ext_exitFlag);
        %c(iType) = c;
        idx = test>=1;
        
        filoInfoFilt = filoInfoBranchStem(idx);
        
        % filter out any that might have passed the exitflag criteria but NOT
        % gave a number for the fit ==0  % maybe flag above later...
        test2 = vertcat(filoInfoFilt(:).Ext_endpointCoordFitPix);
        idx2 = ~isnan(test2);
        filoInfoFilt = filoInfoFilt(idx2);
   value =  nanmean(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic; 
    valueMax = nanmax(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic; 
  text(nx/10,30,['Mean Length ' num2str(value,2) ' um'], 'color',textColor);     
   text(nx/10,50,['Max Length ' num2str(valueMax,2) ' um'],'color',textColor); 
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png']) 
        
        close gcf 
        imageNum = imageNum+1; 
        clear filoInfoFilt
 %% 
h  = setFigure(nx,ny); 
%figure; 

 imshow(-img,[]);   
 hold on 
 
 cellfun(@(x) plot(x(:,2),x(:,1),'y'),edgeYX); 
 filoInfoHigherOrder = filoInfo(type==2); 
 plotfilosIntAndExt(filoInfoHigherOrder,imgSize,1,1,'m',1); 
 
 text(nx/10,10,'High Confidence 1st Order Branches', 'color',textColor);
 pixels = round(10/pixSizeMic); 
   plotScaleBar(pixels,pixels/10,'Label','10um','Color',textColor);
 test = vertcat(filoInfoHigherOrder(:).Ext_exitFlag);
        %c(iType) = c;
        idx = test>=1;
        
        filoInfoFilt = filoInfoHigherOrder(idx);
        
        % filter out any that might have passed the exitflag criteria but NOT
        % gave a number for the fit ==0  % maybe flag above later...
        test2 = vertcat(filoInfoFilt(:).Ext_endpointCoordFitPix);
        idx2 = ~isnan(test2);
        filoInfoFilt = filoInfoFilt(idx2);
   value =  nanmean(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic; 
    valueMax = nanmax(vertcat(filoInfoFilt(:).Ext_length)).*pixSizeMic; 
  text(nx/10,30,['Mean Length ' num2str(value,2) ' um'], 'color',textColor);   
  text(nx/10,50,['Max Length ' num2str(valueMax,2) ' um'],'color',textColor); 
   print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], ...
            [saveDir filesep 'ReconstructMovie' num2str(imageNum) '.png']) 
 
 
    
   
  
        
 
 
 
   cd(saveDir)
 execute = 'mencoder mf://*.png -mf w=800:h=600:fps=0.5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
 system(execute);      
  
 
 
 
 
end 
     
