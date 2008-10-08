function [params,not2keep]=spcklMovCalcBorder(params)
% set the upper, lower, left, and right borders around the images
% to avoid edge effects.  the borders get cropped off when the images are
% saved in spcklMovWrtIm

switch params.nModel
    case 1 %stationary - no substantial border needed
        % add some in case diffusion causes border effect
        params.border.left=10;
        params.border.right=10;
        params.border.top=10;
        params.border.bottom=10;

    case 2 %converging - need border on all sides
        not2keep.nmPerSecFlowSpeed=params.umPerMinFlowSpeed*(1000/60); % convert to nm/s
        nmPerFrame=not2keep.nmPerSecFlowSpeed*params.nSecPerFrame; % calculation flow speed in nm/frame
        b=ceil((nmPerFrame*params.nFrames)/params.pixNM); % max distance covered over movie
        
        params.border.left=b;
        params.border.right=b;
        params.border.top=b;
        params.border.bottom=b;

    case 3 %flow in one direction
        not2keep.nmPerSecFlowSpeed=params.umPerMinFlowSpeed*(1000/60); % convert to nm/s
        nmPerFrame=not2keep.nmPerSecFlowSpeed*params.nSecPerFrame; % calculation flow speed in nm/frame
        b=ceil((nmPerFrame*params.nFrames)/params.pixNM); % max distance covered over movie

        params.border.left=b;
        params.border.right=b;
        params.border.top=b;
        params.border.bottom=b;

    case 4 %MT model - need border on left and right only
        not2keep.nmPerSecFlowSpeed=(params.fastFlowMean+3*params.fastFlowSigma)*(1000/60);
        nmPerFrame=not2keep.nmPerSecFlowSpeed*params.nSecPerFrame; % calculation flow speed in nm/frame
        b=ceil((nmPerFrame*params.nFrames)/params.pixNM); % max distance covered over movie
        
        params.border.left=b;
        params.border.right=b;
        params.border.top=10;
        params.border.bottom=10;

    otherwise
        disp('need to specify border stuff in genFluor if creating new model')
end

% image length and width with border (pixels) 
not2keep.imgLwB = params.imL + params.border.top + params.border.bottom;
not2keep.imgWwB = params.imW + params.border.left + params.border.right;

% image length and width with border (nm) 
not2keep.imgLnm = params.pixNM*not2keep.imgLwB;
not2keep.imgWnm = params.pixNM*not2keep.imgWwB;

% calculate how many sampling squares to use
not2keep.nSamplesL=ceil(not2keep.imgLnm/params.sampleScale);
not2keep.nSamplesW=ceil(not2keep.imgWnm/params.sampleScale);
not2keep.nSampleSquares=not2keep.nSamplesL*not2keep.nSamplesW;
not2keep.sampledIm=zeros(not2keep.nSamplesL,not2keep.nSamplesW);