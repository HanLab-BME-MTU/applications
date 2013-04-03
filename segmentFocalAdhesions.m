function segmentFocalAdhesions(movieData,varargin)


%% --------------- Parameters ---------- %%

pString = 'mask_'; %Prefix for saving masks to file
dName = 'channel ';%Sub-directory name for per-channel masks

%% ------------------ Input ---------------- %%

ip = inputParser;
ip.addRequired('movieData',@(x)(isa(x,'MovieData')));
ip.addOptional('paramsIn',[],@isstruct);

ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('FocalAdhesionSegmentationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(FocalAdhesionSegmentationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end

adSegProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(adSegProc,paramsIn);


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

fStr = ['%0' num2str(floor(log10(nImages))) + 1 '.0f'];%Format string for zero-padding mask file anmes


%Get mask and image directories
maskDirs  = adSegProc.outFilePaths_(p.ChannelIndex);

%TEMP - check for and load thresholding/mask refinement 


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
        

        %Load the current image        
        currImage = movieData.channels_(p.ChannelIndex(iChan)).loadImage(iImage);
        
        
        currMask = blobSegmentThreshold(currImage,0);
        
        
        %write the mask to file                    
        imwrite(currMask,[maskDirs{iChan} filesep pString num2str(iImage,fStr) '.tif']);
        
        if ishandle(wtBar) && mod(iImage,5)
            %Update the waitbar occasionally to minimize slowdown
            waitbar((iImage + (iChan-1)*nImages) / nImTot,wtBar)
        end
                
    
    end    
   
    
end

if ishandle(wtBar)
    close(wtBar)
end

%% ---------------- Finalization, Output --------------- %%

adSegProc.setDateTime;
movieData.save; %Save the new movieData to disk


