function [analInfo,saveTimeVect ] = GCAReconstructFilopodia(imDir,framesPath,troubleshoot,sortOn,restart,analInfo,startFrameIn,protrusion)


 
%% GCAReconstructFilopodia - changed from reconstructFilopodia 20141022
% This function follows neurite body estimation and applies a small scale 
% steerable filter to detect filopodia ridges. A graph matching based
% method is then used to resconstruct the entire neurite body structure including 
% .
%
% detection params.: 
%                 .filterParams.filterOrder
%                 .filterParams.scales
% These were previous parameters from then neurite body estimation (NOTE
% 06-09 trying to separate these into separate steps and functions) 
% workflow will be such that 1) smooth body estimation, 2) protrusion
% vector calc ,3) reconstruction (current function) including orientation
% filopodia body orientation calculations using the protrusionVectors
%                 .localThresh = 1 if you want to use local thresholding
%                 .saveViaBB = 1 if you want to save via backbone 
%                 .localThresh.patchSize = patch size for thresholding 
%
% OUTPUT: adds a field to analInfo called filoInfo
% filoInfo is a N structure x 1 structure with fields providing information
% regarding each filopodia where N is the number of filopodia pieces
% reconstructed in iFrame.  Certain filopodia are clustered into branch
% groups via the field groupCount (though there might be more elegant ways
% to store this information) Nested fields were avoided here to facilitate the data
% extraction in later steps (ie each filopodia was given an field ID specifying branch order rather than nesting
% the structure)

%% Check Input% NOTE TO SELF FIX INPUT 
filterOrder = 4; 
scales = 1.5; 

  %filterOrder = detectionParams.filterParams.filterOrder;  
  %scales = detectionParams.filterParams.scales; 
  %save([framesPath filesep 'DetectionParams.mat'],'detectionParams'); 

%save([framesPath filesep 'detectionParams.mat'],'detectionParams'); 

% load protrusionVectors from the body
if ~isempty(protrusion)
% load([protrusion filesep 'protrusion_vectors.mat']); 
 normals = protrusion.normals; % need to load the normal input from the smoothed edges in order  to calculation the filopodia body orientation
 smoothedEdge = protrusion.smoothedEdge; 
end

%% load movie and prepare output movie options 

if restart ==1;
  %numFinished = numel(analInfo);
%    
% 
 
startFrame = find(arrayfun(@(x) isempty(x.filoInfo),analInfo),1,'first')-1; 
display(['Restarting at Frame ' num2str(startFrame)]); 
% startFrame = startFrame-1;
%startFrame = 81;
  %  startFrame = numFinished; 
  elseif restart ==2
      startFrame = startFrameIn; 

else startFrame = 1; 
end 


% collect images and initiate
[listOfImages] = searchFiles('.tif',[],imDir,0);
nImTot = size(listOfImages,1);

if sortOn ==1 % reminder make it check for padding automatically! 

 for iFrame =  1:nImTot-1;
        imageName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
        
        %Call "getFilenameBody" from common dir to split filename into path,
        %body, no, and ext. Put path, body and ext into appropriate cell/number vector for that
        %frame
        [path body no ext ] = getFilenameBody(imageName);
        
        
        pathCell(iFrame) = cellstr(path);
        bodyCell(iFrame) = cellstr(body);
        extCell(iFrame) = cellstr(ext);
        
        % output "no" is a string so convert to number to sort
        num(iFrame)  = str2double(no);
        
  end
    %Sort number vector numerically
    sortednum = sort(num);
    
    %Convert number vector to cell
    sortednum_cell = num2cell(sortednum);
    
    %Create Sorted Image List
    sortedImages = [pathCell', bodyCell', sortednum_cell', extCell'];
%     if restart ==1; 
%         sortedImages = sortedImages(startFrame:end,:);
%     end 
end 
   
 
%% Configure Figure 
if sortOn~=1        
 fileName = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
else         
          fileName = [char(sortedImages(1,1)) filesep char(sortedImages(1,2)),...;
                  num2str(sortednum(1)) char(sortedImages(1,4))];
end 

              img = double(imread(fileName)); 
        [ny nx] = size(img) ; 
        zoom = 1; 
        h = setFigure(nx,ny,'off'); 
        
        fmt = ['%0' num2str(ceil(log10(nImTot))) 'd'];

        framesPath = [framesPath filesep]; 

        ext = '.png';

        [numFrames ~ ] = size(listOfImages);
        
%% Start Loop         
for iFrame = startFrame:nImTot-1
  %try
        %diary([framesPath filesep 'Bug_Frame' num2str(iFrame) '.txt'])
        
   %iFrame = 1; % just do for one frame for now
if sortOn~=1
    fileName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
else 
    fileName = [char(sortedImages(iFrame,1)) filesep char(sortedImages(iFrame,2)),...;
        num2str(sortednum(iFrame)) char(sortedImages(iFrame,4))];
end   
analInfo(iFrame).imgPointer = fileName; % right now this is just useful for me 
% should REALLY try to familialize myself with movieData!! 
   img = double(imread(fileName)); 
 
    %% reset fig (shouldn't have to do this but having some trouble resetting the handle)
            h = setFigure(nx,ny,'off'); 
%% Background Fluorescence Filter (3*first mode)
%perform an initial thresholding first to get out main components 
[maskBack,backMu,backSig] = estimateBackgroundArea(img); 
% imgToThresh = img.*~maskBack; % note doing this is not so awesome for 
% for the local thresholding as if you input this the local threshold will
% then begin picking up on these gradients  Note 12-15 changed

%% Multi scale ridge detect small scales
    
[maxRes, maxTh ,maxNMS ,scaleMap]= multiscaleSteerableDetector(img,filterOrder,scales); 

% save these as we will need them as input for later 
analInfo(iFrame).filterInfo.maxTh = maxTh;
analInfo(iFrame).filterInfo.maxRes = maxRes;

analInfo(iFrame).scaleMap = scaleMap; 

%% New thresholding of small scale ridge NMS
 forValues = maxNMS.*~maskBack; % take out background response based on fluorescence intensity
 valuesFilter = forValues(forValues~=0); 
 

    % try a second step where you filter out the background response from the 
    % estimation around the area of interest: keep thresholding on a pixel
    % by pixel basis. 
    
    [respNMSMean,respNMSSTD]   = fitGaussianModeToPDF(valuesFilter); 
    cutoffTrueResponse = respNMSMean+3*respNMSSTD; 
    n1 = hist(valuesFilter,100);
    
        if troubleshoot ==1 % plot the histogram with cut-off overlay so can see what losing 
        ridgeDiagDir = [framesPath filesep 'candidateRidges'];
        if ~isdir(ridgeDiagDir)
         mkdir(ridgeDiagDir)
        end 

        hist(valuesFilter,100); 
        hold on ed
        line([cutoffTrueResponse cutoffTrueResponse],[0,max(n1)],'color','r','Linewidth',2); 
        title('Red line 3*std of first mode'); 
        saveas(gcf,[ridgeDiagDir filesep 'maxNMSResFirstModeCut' num2str(iFrame,fmt) '.eps'],'psc2');
        close gcf
        end
        
    canRidges = maxNMS.*~maskBack;
    canRidges(canRidges<cutoffTrueResponse) = 0; 
    
    analInfo(iFrame).filterInfo.ThreshNMS = canRidges; 
  

%% load neurite body estimation
erodfilo = analInfo(iFrame).masks.neuriteEdge ;  
         
%% clean the thresholded NMS response using the body estimation information 
    
skelIn = canRidges; 

skelIn = bwmorph(skelIn,'thin',inf ); 


% Clean up ridges (ie potential filo) by looking getting ride of those that
% are within the "neurite body"
% Here defined by erosion of locally thresholded image. 
input = analInfo(iFrame) ; % NOTE Eventually just make the input a structure easier to deal with. 

if ~isempty(protrusion) 
normalC = normals{iFrame}; 
smoothedEdgeC = smoothedEdge{iFrame}; % FIX THIS BACK 
else  
    normalC = []; 
    smoothedEdgeC = []; 
    
end 
%cd(framesPath);
% Main Function (should rename) that performs the reconstructions
[reconstruct,filoInfo] = gcaAttachFilopodiaStructures(img,skelIn,erodfilo,scaleMap,maxRes,maxTh,nImTot,input,normalC,smoothedEdgeC,iFrame,framesPath);

analInfo(iFrame).filoInfo = filoInfo; % will already be filtered for short filo

analInfo(iFrame).reconstructInfo = reconstruct;

analInfo(iFrame).reconstructInfo.createTime = clock; 

tic 
save([framesPath filesep 'analInfoTestSave.mat'],'analInfo','-v7.3'); 
toc
saveTimeVect(iFrame) = toc; 

makePlots =1;

%% Plot the Fits  
if makePlots == 1
     filoCoordsDir = [framesPath filesep 'Final Fits']; 
 subDir{1} = [filoCoordsDir filesep 'WithInternal']; 
 subDir{2} = [filoCoordsDir filesep 'NoInternal']; 
 
 if ~isdir([filoCoordsDir filesep 'WithInternal']) 
     mkdir([filoCoordsDir filesep 'WithInternal'])
 end 
 
 if ~isdir([filoCoordsDir filesep 'NoInternal']) 
     mkdir([filoCoordsDir filesep 'NoInternal'])
 end 
 img = [img img]; 
 for iDir = 1:2
     
     [nyD,nxD] = size(img); 
     h = setFigure(nxD,nyD,'off');
     imshow(-img,[]) ;
     hold on
     
     neuriteEdge = bwboundaries(erodfilo);
     cellfun(@(x) plot(x(:,2),x(:,1),'y'),neuriteEdge);
     hold on
     if iDir ==1;
         justExtFlag = 0;
         plotOrientFlag = 0;% donr't plot orient for internal 
     else
         justExtFlag = 1;
         plotOrientFlag =1; % 
     end
    % plotOrientFlag = 1; 
     plotfilosIntAndExt(filoInfo,[ny,nx],1,justExtFlag,[],plotOrientFlag); % Note small bug here if don't filter 
     % by fit Dimensions of matrices being concatenated are not consistent.
     % line 39 try to fix later!! 2013_07_14
     secPerFrame = 5; 
     pixSizeMic = 0.216; 
     
     text(10,10,[num2str(iFrame*secPerFrame-secPerFrame) ' s'],'Color','k');
     pixels = round(10/pixSizeMic); 
     
     plotScaleBar(pixels,pixels/10,'Label','10 um','Color',[0 0 0]); 
     ext = '.png';
     print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], [subDir{iDir} filesep 'finalFits' num2str(iFrame, fmt) ext])
     saveas(h,[subDir{iDir} filesep 'finalFits' num2str(iFrame,fmt) '.eps'],'psc2');
     close gcf
 end % iDir
end % if makeplots
    
end %end iframe 

for iDir = 1:2
 cd(subDir{iDir})
 execute = 'mencoder mf://*.png -mf w=800:h=600:fps=5:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o movie.wmv';
 system(execute); 
end 



clear all

end % 


        
        
        




