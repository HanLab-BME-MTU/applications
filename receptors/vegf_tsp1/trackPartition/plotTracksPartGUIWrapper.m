function plotTracksPartGUIWrapper(trackMD,options)

% Get channels
process = trackMD.processes_{trackMD.getProcessIndex('TrackPartitionProcess')};
trackChannel = process.trackChannel_;
maskChannel = process.maskChannel_;
maskMDPath = process.maskMovieDataPath_;
maskMDFileName = process.maskMovieDataFileName_;

% Load mask MD
maskMD = load([maskMDPath,filesep,maskMDFileName]);
maskMD = maskMD.MD;

% Get partitioned tracks
output = load(process.outFilePaths_{trackChannel});
tracks = output.tracksPart;
diffAnalysisRes = output.diffAnalysisRes;

if options.showImage
    % Load the image
    trackChannelInfo = trackMD.channels_(trackChannel);
    trackImage = imread(trackChannelInfo.channelPath_);
    maskChannelInfo = maskMD.channels_(maskChannel);
    maskImage = imread(maskChannelInfo.channelPath_);
    trackColor = parseColor(trackChannelInfo);
    maskColor = parseColor(maskChannelInfo);
    
    trackImage = mat2gray(trackImage);
    maskImage = mat2gray(maskImage);
    
    if options.whiteImage 
        % Create RGB image subtractively
        I = ones([size(trackImage),3]);
        switch trackColor
            case 'blue'
                I(:,:,1) = I(:,:,1)-trackImage;
                I(:,:,2) = I(:,:,2)-trackImage;
            case 'green'
                I(:,:,1) = I(:,:,1)-trackImage;
                I(:,:,3) = I(:,:,3)-trackImage;
            case 'red'
                I(:,:,2) = I(:,:,2)-trackImage;
                I(:,:,3) = I(:,:,3)-trackImage;
        end
        switch maskColor
            case 'blue'
                I(:,:,1) = I(:,:,1)-maskImage;
                I(:,:,2) = I(:,:,2)-maskImage;
            case 'green'
                I(:,:,1) = I(:,:,1)-maskImage;
                I(:,:,3) = I(:,:,3)-maskImage;
            case 'red'
                I(:,:,2) = I(:,:,2)-maskImage;
                I(:,:,3) = I(:,:,3)-maskImage;
        end
    else
        % Create RGB image additively
        I = zeros([size(trackImage),3]);
        switch trackColor
            case 'blue'
                I(:,:,3) = trackImage;
            case 'green'
                I(:,:,2) = trackImage;
            case 'red'
                I(:,:,1) = trackImage;
        end
        switch maskColor
            case 'blue'
                I(:,:,3) = maskImage;
            case 'green'
                I(:,:,2) = maskImage;
            case 'red'
                I(:,:,1) = maskImage;
        end
    end
else
    I = [];
end

timeRange = options.timeRange;
newFigure = 1;
if options.showConf 
    showConf = 1;
else 
    showConf = 0;
end
simplifyLin = options.simplifyLin;
offset = [];
plotSubset = options.plotSubset;
if options.showBoxes
    showBoxes = 1;
else
    showBoxes = 0;
end

plotTracksPart(tracks,diffAnalysisRes,timeRange,newFigure,I,showConf,...
    simplifyLin,offset,plotSubset,showBoxes);
    


end

function color = parseColor(channelInfo)
% Return a primary color based on channel wavelength
excitation = channelInfo.excitationWavelength_;
emission = channelInfo.emissionWavelength_;
if ~isempty(excitation)
    if excitation < 460
        color = 'blue';
    elseif (excitation >= 460) && (excitation < 550)
        color = 'green';
    else excitation >= 550
        color = 'red';
    end
else
    if emission < 490
        color = 'blue';
    elseif (emission >= 490) && (emission < 580)
        color = 'green';
    else 
        color = 'red';
    end
end
end