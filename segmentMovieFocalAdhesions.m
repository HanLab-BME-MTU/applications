function segmentMovieFocalAdhesions(movieData)
%SEGMENTFOCALADHESIONS segments and labels focal adhesions in the input movie
%
% segmentFocalAdhesions(movieData)
%
% This function uses blobSegmentThreshold to get an initial segmentation of
% adhesions, and then refines this segmentation in two ways:
%   
%   Splitting of large groups of focal adhesions using steerable filtering:
%       Closely spaced clusters of adhesions can be merged by blobSegment,
%       so this splits these clusters into individual adhesions by
%       identifying troughs in the intensity via steerable filtering.
%
%   Morphological post-processing and labelling:
%       Removes adhesions which are either very small or persist for a very
%       short duration by specifying a minimum area*volume for a  segmented
%       adhesion. 
%       Also minimizes spurious merging of closely apposed adhesions by
%       performing a morphological opening in time and space, to remove
%       connections between adhesions which are either very thin, very
%       brief, or both.
%       Finally, each adhesions which is contiguous in space & time (e.g.
%       overlaps in each consecutive frame) is assigned a unique integer
%       label.
%


%Hunter Elliott
%April 2013


%% --------------- Parameters ---------- %%

pString = 'mask_'; %Prefix for saving masks to file
dName = 'channel ';%Sub-directory name for per-channel masks

sfOrd = 2;%Order of steerable filters to use for splitting adjacent adhesions. MUST BE EVEN ORDER (2 or 4)
connNum = 6;%Connectivity to use for labelling 3D mask. We use lowest to minimize unintentional merging.

%% ------------------ Input ---------------- %%

ip = inputParser;
ip.addRequired('movieData',@(x)(isa(x,'MovieData')));
ip.parse(movieData);

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('FocalAdhesionSegmentationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    error('No FocalAdhesionSegmentationProcess in input movieData! please create the process and use the process.run method to run this function!')
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(FocalAdhesionSegmentationProcess(movieData,...
%         movieData.outputDirectory_));                                                                                                 
end

adSegProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(adSegProc);


%% --------------- Init -------------- %%

nChanSeg = numel(p.ChannelIndex);


% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex    
    inFilePaths{1,i} = movieData.getChannelPaths{i};    
end
adSegProc.setInFilePaths(inFilePaths);

%Set up the mask directories as sub-directories of the output directory
for j = 1:nChanSeg;
    
    %Create string for current directory
    currDir = [p.OutputDirectory filesep dName num2str(p.ChannelIndex(j))];    
    %Save this in the process object
    adSegProc.setOutMaskPath(p.ChannelIndex(j),currDir);
   
    %Check/create directory
    mkClrDir(currDir)               
end



nImages = movieData.nFrames_;   
nImTot = nImages * nChanSeg;
imSize = movieData.imSize_;

fStr = ['%0' num2str(floor(log10(nImages))) + 1 '.0f'];%Format string for zero-padding mask file anmes


%Get mask and image directories
maskDirs  = adSegProc.outFilePaths_(p.ChannelIndex);

%TEMP - check for and load thresholding/mask refinement 

masks = false([imSize, nImages nChanSeg]);

%% ------------ Segmentation -------------- %%

if ~p.BatchMode && feature('ShowFigureWindows')
    wtBar = waitbar(0,['Please wait, segmenting adhesions in channel ' num2str(p.ChannelIndex(1)) ' ...']);        
else
    wtBar = -1;
end        


for iChan = 1:nChanSeg
        
        
    if ishandle(wtBar)        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, segmenting adhesions in channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end    
    
    
    for iImage = 1:nImages
        

        % ---------- Initial Segmentation ------------- %
        
        %Load the current image        
        currImage = movieData.channels_(p.ChannelIndex(iChan)).loadImage(iImage);
        
        %Run the blob segmentation to get approx outline. This removes
        %background well but tends to merge close adhesions        
        currMask = blobSegmentThreshold(currImage,0);                

        
        % --------- Splitting of Adjacent Adhesions --------- %
        %Uses steerable filtering to split larger adhesions based on
        %troughs in intensity within them.
        
        %If selected, split adjacent adhesions using steerable filtering
        if p.SteerableFilterSigma > 0
            
            %Convert to pixels, with minimum of 1.
            filtSigma = max(p.SteerableFilterSigma / movieData.pixelSize_,1);
           
            %Filter with steerable filter that responds to ridge-like
            %objects
            imFilt = steerableDetector(double(currImage),sfOrd,filtSigma);
            
            %Since the filter will also respond to troughs in intensity by
            %orienting perpendicular to them, we compare this to the
            %response to the inverted image. 
            imInvFilt = steerableDetector(double(max(currImage(:)) - currImage),sfOrd,filtSigma);
            
            %Regions where the response is stronger in the inverted image
            %are troughs and are then removed.
            currMask(imInvFilt > imFilt) = false;            
            
        end
        
        %Store this in the combined matrix at the end, to make indexing
        %operations above simpler
        masks(:,:,iImage,iChan) = currMask;
        
         
        if ishandle(wtBar) && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end
                
    
    end    
   
    
end



%% ------------ Morphological Post-Processing ------------ %%
%Splits adjacent adhesions which are weakly connected using mathematical
%morphology

if ishandle(wtBar)        
    waitbar(0,wtBar,['Please wait, post-processing segmented adhesions in channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
end    
    
%HLE - this can be re-written to do all channels simultaneously. Might be
%slightly faster, but might use more memory as well...

for iChan = 1:nChanSeg    
    
    if p.MinVolTime > 0
        
        %Convert to pixel^2*frames
        minSize = ceil(p.MinVolTime / ((movieData.pixelSize_/1e3)^2 * movieData.timeInterval_));
        
        %Remove the small and / or short lived objects
        masks(:,:,:,iChan) = bwareaopen(masks(:,:,:,iChan),minSize,connNum);        
        
    end
    
        
    if p.OpeningRadiusXY > 0 || p.OpeningHeightT > 0

        %Convert radii into pixels & frames
        rOpenXY = ceil(p.OpeningRadiusXY / movieData.pixelSize_);
        rOpenT = ceil(p.OpeningHeightT / movieData.timeInterval_);

        nHoodOpen = strel('disk',rOpenXY,0);
        nHoodOpen = nHoodOpen.getnhood;

        if rOpenT > 0
            nHoodOpen = repmat(nHoodOpen,[1 1 rOpenT]);
        end

        %Perform the opening on the 3D mask matrix. The spatial (XY) radius will
        %split adhesions which are connected only by very thin strips of
        %pixels, while the height in time will remove connections which are
        %thicker but that only persist for a short period of time.
        masks(:,:,:,iChan) = imopen(masks(:,:,:,iChan),nHoodOpen);
       
    end                    

    
    if ishandle(wtBar)        
        waitbar(iChan/nChanSeg,wtBar)
    end
    
end

%Label each contiguous object. This is the "tracking" step, and also allows
%us to visualize which adhesions are merged. We do it the two-step way to
%minimize memory usage with these potentially large matrices. We also do
%this on the full multi-channel matrix so that objects are labelled
%uniquely between the channels (by specifying a 3D connectivity, objects in
%diff channels are not considered adjacent)
cc = bwconncomp(masks,connNum);
masks = labelmatrix(cc);



%% ------------- Output ----------------- %%
%Writes the final masks to disk



for iChan = 1:nChanSeg
    
    if ishandle(wtBar)        
        waitbar((iChan-1)*nImages / nImTot,wtBar,['Please wait, saving masks to disk for channel ' num2str(p.ChannelIndex(iChan)) ' ...']);        
    end   
    
    for iImage = 1:nImages

        %write the mask to file                    
        imwrite(masks(:,:,iImage,iChan),[maskDirs{iChan} filesep pString num2str(iImage,fStr) '.tif']);
        
        if ishandle(wtBar) && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end

    end
end

%Store the maximum number of objects in the process class
adSegProc.setMaxIndex(cc.NumObjects);


if ishandle(wtBar)
    close(wtBar)
end




