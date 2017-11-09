function [out,jump] = fixWindowJumpsNew(movieObj,varargin)
%This function tries to corrects window's jumps
%Usage:
%
%      [out,jump] = fixWindowJumps(movieObj,inMap)
%
%Input:
%
%      movieObj         - movie list or movie data object
%      inMap(optional)  - additional map to be corrected (e.g. biosensor activity)
%      samplingProcess  - either 'WindowSamplingProcess'(default) or 'ProtrusionSamplingProcess'
%Output:
%
%       out - corrected map
%       jump - structure with the location and magnitude of the corrected
%       jumps
%Marco Vilela, 2014

ip = inputParser;
ip.addRequired('movieObj',@(x) isa(x,'MovieList') || isa(x,'MovieData'));


if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

ip.addParameter('inMap',[]);
ip.addParameter('samplingProcess','WindowSamplingProcess');
ip.addParameter('signalChannel',1,@isscalar);

ip.parse(movieObj,varargin{:});
inMap            = ip.Results.inMap;
samplingProcess  = ip.Results.samplingProcess;
channel          = ip.Results.signalChannel;


nCell      = numel(ML.movies_);
out        = cell(1,nCell);
jump       = cell(1,nCell);
minPixMove = 8;

for iCell = 1:nCell
    
    currMD = ML.movies_{iCell};
    idx    = currMD.getProcessIndex(samplingProcess);
    
    if strcmp(samplingProcess,'WindowSamplingProcess')
        map    = squeeze(currMD.processes_{idx}.loadChannelOutput(channel).avg(:,1,:));
    else
        map    = currMD.processes_{idx}.loadChannelOutput.avgNormal;
    end
    
    if isempty(inMap)
        mapAux = map;
    else
        mapAux = inMap;
    end
    
    out{iCell} = mapAux;
    distP      = nan(1,currMD.nFrames_);
    idx        = currMD.getProcessIndex('WindowingProcess');
    
    for iF = 1:currMD.nFrames_
        startW = currMD.processes_{idx}.loadChannelOutput(iF);
        iC     = 1;
        if isempty(startW{iC})
            while isempty(startW{iC})
                iC = iC + 1;
            end
        else
            while isempty(startW{iC}{1})
                iC = iC + 1;
            end
        end
            
        startP    = startW{iC}{1}{end}(:,end);
        distP(iF) = norm(startP);
    end
    shiftP        = diff(distP);
    
    if any(shiftP > minPixMove)
        
        breaks         = detectOutliers(shiftP,minPixMove);
        [outMap,fJump] = forwardShift(breaks,shiftP,mapAux,currMD);
        jump{iCell}    = fJump;
        out{iCell}     = outMap;
        
    else
        
        out{iCell}     = mapAux;
        jump{iCell}.location  = [];
        jump{iCell}.magnitude = [];
        
    end
    
    
end
clear distP
clear xCorr

end


function xCorr = circularNanCrossCor(x,y)

nPoint  = numel(y);
circIdx = gallery('circul',nPoint);

xCorr = nan(1,nPoint);

if ~(sum(isnan(y)) == nPoint || sum(isnan(x)) == nPoint )
    
    for i = 1:nPoint
        
        xCorr(i) = nanCrossCorrelation( x, y(circIdx(i,:)) );
        
    end
    
end

end

function [outMap,jump] = forwardShift(breaks,shiftP,mapAux,currMD)

jump.location  = [];
jump.magnitude = [];
outMap         = mapAux;

if ~isempty(breaks)
    bPeaks     = findBlock(breaks,1);
    nB         = numel(bPeaks);
    singleJump = [];
    blockJump  = [];
    map        = mapAux;
    
    for iB = 1:nB
        if length(bPeaks{iB}) == 2
            singleJump = [singleJump;bPeaks{iB}(1:2:numel(bPeaks{iB}))+1];
        else
            blockJump = [blockJump;bPeaks{iB}];
        end
    end
    
    
    jC = 1;
    
    if ~isempty(singleJump)
        for iJ = 1:numel(singleJump)
            ccIdx = 1;
            while sum(isfinite(map(:,singleJump(iJ)-ccIdx))) < 40%Minimum number of points for xCorr
                ccIdx = ccIdx + 1;
                if singleJump(iJ)-ccIdx < 1
                    break;
                end
            end
            
            xCorr    = circularNanCrossCor(outMap(:,singleJump(iJ)),outMap(:,singleJump(iJ)-ccIdx));
            [~,idx1] = max(xCorr);
            
            outMap(:,singleJump(iJ)) = circshift(outMap(:,singleJump(iJ)),-idx1);
            jump.location{jC}        = singleJump(iJ);
            jump.magnitude{jC}       = -idx1;
            jC                       = jC + 1;
            
            if abs( abs(shiftP(singleJump(iJ))) - abs(shiftP(singleJump(iJ)-1)) ) > 4
                
                xCorr    = circularNanCrossCor(outMap(:,singleJump(iJ)+1),outMap(:,singleJump(iJ)));
                [~,idx1] = max(xCorr);
                
                outMap(:,singleJump(iJ)+1:end) = circshift(mapAux(:,singleJump(iJ)+1:end),-idx1);
                jump.location{jC}        = singleJump(iJ)+1:currMD.nFrames_;
                jump.magnitude{jC}       = -idx1;
                jC                       = jC + 1;
                
            end
            
        end
        
    end
    
    if ~isempty(blockJump)
        
        for iJ = 1:numel(blockJump) - 1
            
            moveW = blockJump(iJ)+1:currMD.nFrames_;
            ccIdx = 1;
            while sum(isfinite(outMap(:,moveW(1)-ccIdx))) < 40%Minimum number of points for xCorr
                ccIdx = ccIdx + 1;
                if ccIdx < 0
                    break;
                end
            end
            xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)-ccIdx));
            [~,idx1] = max(xCorr);
            
            outMap(:,moveW)    = circshift(outMap(:,moveW),-idx1);
            jump.location{jC}  = moveW;
            jump.magnitude{jC} = -idx1;
            jC                 = jC + 1;
            
        end
        moveW = blockJump(end)+1:currMD.nFrames_;
        ccIdx = 1;
        while sum(isfinite(outMap(:,moveW(1)-ccIdx))) < 40%Minimum number of points for xCorr
            ccIdx = ccIdx + 1;
            if ccIdx < 0
                break;
            end
        end
        xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)-ccIdx));
        [~,idx1] = max(xCorr);
        
        outMap(:,moveW)    = circshift(outMap(:,moveW),-idx1);
        jump.location{jC}  = moveW;
        jump.magnitude{jC} = -idx1;
        jC                 = jC + 1;
    end
    
    
end%

end