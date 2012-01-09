function analyze3DMovieMaskedIntensities(movieData,varargin)

%TEMP OOOOOOOBBBBBBBVVVVVIIIOOOUUUSSSSSLLLLYYYYY
p.ChannelIndex = 1:3;
p.BatchMode = false;

%Make sure the movie has been segmented and has masks
iSegProc = movieData.getProcessIndex('SegmentationProcess3D',1,~p.BatchMode);

if isempty(iSegProc) || ~movieData.processes_{iSegProc}.checkChannelOutput(p.ChannelIndex(1))%TEMP? how do deal with three-input channels for mask and one output?? ALways associate with channel 1?
    error('The input MovieData does not have valid masks for the selected channel!');
end

%Get mask file names and directory
maskDir = movieData.processes_{iSegProc}.outFilePaths_{1,p.ChannelIndex(1)};
maskNames = movieData.processes_{iSegProc}.getOutMaskFileNames(p.ChannelIndex(1));


imDir = movieData.getChannelPaths(p.ChannelIndex);
imNames = movieData.getImageFileNames(p.ChannelIndex);

nChan = numel(p.ChannelIndex);
nFrames = movieData.nFrames_;

imSize = [movieData.imSize_ movieData.nSlices_];

for iFrame = 1:nFrames
    
    
    for iChan = 1:nChan    
        if iChan == 1
            currIm = zeros([imSize nChan],'uint16');
        end
        currIm(:,:,:,iChan) = stackRead([imDir{iChan} filesep imNames{iChan}{iFrame}]);                                    
    end
    
    currMask = tif3Dread([maskDir filesep maskNames{1}{iFrame}]);
    
    branchProfiles = analyze3DImageMaskedIntensities(currIm,currMask);
    
    
end
