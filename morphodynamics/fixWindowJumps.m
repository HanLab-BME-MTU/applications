function [out,jump] = fixWindowJumps(movieObj,varargin)
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


nCell = numel(ML.movies_);
out   = cell(1,nCell);
jump  = cell(1,nCell);

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
        while isempty(startW{iC})
            iC = iC + 1;
        end
        startP    = startW{iC}{1}{end}(:,end);
        distP(iF) = norm(startP);
    end
    shiftP        = diff(distP);
    
    if any(shiftP > 8)
        breaks         = detectOutliers(shiftP,8);
        [bb,bl]        = findBlock(setdiff(1:currMD.nFrames_,breaks),1);
        [~,idxBl]      = max(bl);
        
        breaksForw     = breaks(breaks > bb{idxBl}(1));
        
        [outMap,fJump] = forwardShift(breaksForw,mapAux,currMD,shiftP(breaksForw));
        revMap         = fliplr(outMap);
        
        breaksBack     = sort(currMD.nFrames_ - ( breaks(breaks < bb{idxBl}(1))));
        [outMap,bJump] = forwardShift(breaksBack,revMap,currMD,shiftP(breaks(breaks < bb{idxBl}(1))));
        if ~isempty(bJump.location)
            bJump.location = cellfun(@(x) 1-x+currMD.nFrames_,bJump.location,'Unif',0);
        end
        jump{iCell}    = catstruct(bJump,fJump);%MergeStruct(bJump,fJump);
        out{iCell}     = fliplr(outMap);
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
for i = 1:nPoint
    
    xCorr(i) = nanCrossCorrelation( x, y(circIdx(i,:)) );
    
end

end

function [outMap,jump] = forwardShift(breaks,mapAux,currMD,shiftP)

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
                if ccIdx < 0
                    break;
                end
            end
            
            xCorr    = circularNanCrossCor(map(:,singleJump(iJ)),map(:,singleJump(iJ)-ccIdx));
            [~,idx1] = max(xCorr);
            
            outMap(:,singleJump(iJ)) = circshift(mapAux(:,singleJump(iJ)),-idx1);
            jump.location{jC}        = singleJump(iJ);
            jump.magnitude{jC}       = -idx1;
            jC                       = jC + 1;
        end
        
        moveW    = singleJump(end)+1:currMD.nFrames_;
        ccIdx = 1;
        while sum(isfinite(outMap(:,moveW(1)+ccIdx))) < 40%Minimum number of points for xCorr
            ccIdx = ccIdx + 1;
            if ccIdx < 0
                break;
            end
        end
        
        xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)+ccIdx));
        [~,idx1] = max(xCorr);
        
        outMap(:,moveW)    = circshift(outMap(:,moveW),-idx1);
        jump.location{jC}  = moveW;
        jump.magnitude{jC} = -idx1;
        jC                 = jC + 1;
    end
    
    if ~isempty(blockJump)
        
        if sum(shiftP < 0) == sum(shiftP > 0)
            
            for iJ = 1:numel(blockJump)/2
                
                moveW    = blockJump(2*iJ-1)+1:blockJump(iJ*2);
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
        else
            
            testSign = sign(shiftP);
            if numel(blockJump) == 1
                
                moveW    = blockJump(1)+1:currMD.nFrames_;
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
                
                
            elseif ( abs(sum(testSign(1:2:end-2))) == abs(sum(testSign(2:2:end-1))) && abs(sum(testSign(2:2:end-1))) == 0.5*(numel(testSign)-1) )
                
                for iJ = 1:floor(numel(blockJump)/2)
                    
                    moveW    = blockJump(2*iJ-1)+1:blockJump(iJ*2);
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
                
                moveW    = blockJump(end)+1:currMD.nFrames_;
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
                
            else
                
                for iJ = blockJump(1):blockJump(end)%floor(numel(blockJump)/2)
                    ccIdx = 1;
                    while sum(isfinite(outMap(:,iJ+ccIdx))) < 40%Minimum number of points for xCorr
                        ccIdx = ccIdx + 1;
                        if ccIdx < 0
                            break;
                        end
                    end
                    xCorr    = circularNanCrossCor(outMap(:,iJ+ccIdx),outMap(:,iJ));
                    [~,idx1] = max(xCorr);
                    
                    outMap(:,iJ+1)     = circshift(outMap(:,iJ+ccIdx),-idx1+1);
                    jump.location{jC}  = iJ;
                    jump.magnitude{jC} = -idx1+1;
                    jC                 = jC + 1;
                    
                end
                
            end
            
            
            
            moveW    = blockJump(end)+1:currMD.nFrames_;
            ccIdx = 1;
            while sum(isfinite(outMap(:,moveW(1)+ccIdx))) < 40%Minimum number of points for xCorr
                ccIdx = ccIdx + 1;
                if ccIdx < 0
                    break;
                end
            end
            xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)+ccIdx));
            [~,idx1] = max(xCorr);
            
            outMap(:,moveW)    = circshift(outMap(:,moveW),idx1);
            jump.location{jC}  = moveW;
            jump.magnitude{jC} = idx1;
            jC                 = jC + 1;
            
            
            
        end
    end
    
end

end