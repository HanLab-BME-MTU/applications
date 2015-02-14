function [out,jump] = fixWindowJumps(movieObj,inMap)
%This function tries to corrects window's jumps  
%Usage:
%
%      [out,jump] = fixWindowJumps(movieObj,inMap)
%
%Input:
%
%      movieObj         - movie list or movie data object
%      inMap(optional)  - additional map to be corrected (e.g. biosensor activity)
%Output:
%
%       out - corrected map
%       jump - structure with the location and magnitude of the corrected
%       jumps
%Marco Vilela, 2014



if isa(movieObj,'MovieData')
    
    ML = movieData2movieList(movieObj);
    
else
    
    ML = movieObj;
    
end

nCell = numel(ML.movies_);
out   = cell(1,nCell);
jump  = cell(1,nCell);

for iCell = 1:nCell
    
    currMD = ML.movies_{iCell};
    idx    = currMD.getProcessIndex('WindowSamplingProcess');
    map    = squeeze(currMD.processes_{idx}.loadChannelOutput(1).avg(:,1,:));
    
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
        jump{iCell}    = MergeStruct(bJump,fJump);
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
    for iJ = 1:numel(singleJump)
        xCorr    = circularNanCrossCor(map(:,singleJump(iJ)),map(:,singleJump(iJ)-1));
        [~,idx1] = max(xCorr);
        
        outMap(:,singleJump(iJ)) = circshift(mapAux(:,singleJump(iJ)),-idx1);
        jump.location{jC}        = singleJump(iJ);
        jump.magnitude{jC}       = -idx1;
        jC                       = jC + 1;
    end
    
    if ~isempty(blockJump)
        
        if sum(shiftP < 0) == sum(shiftP > 0)
            
            for iJ = 1:numel(blockJump)/2
                
                moveW    = blockJump(2*iJ-1)+1:blockJump(iJ*2);
                xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)-1));
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
                xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)-1));
                [~,idx1] = max(xCorr);
                
                outMap(:,moveW)    = circshift(outMap(:,moveW),-idx1);
                jump.location{jC}  = moveW;
                jump.magnitude{jC} = -idx1;
                jC                 = jC + 1;
                
            elseif ( abs(sum(testSign(1:2:end-2))) == abs(sum(testSign(2:2:end-1))) && abs(sum(testSign(2:2:end-1))) == 0.5*(numel(testSign)-1) )
                
                for iJ = 1:floor(numel(blockJump)/2)
                    
                   moveW    = blockJump(2*iJ-1)+1:blockJump(iJ*2);
                   xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)-1));
                   [~,idx1] = max(abs(xCorr));
                   
                   outMap(:,moveW)    = circshift(outMap(:,moveW),-idx1);
                   jump.location{jC}  = moveW;
                   jump.magnitude{jC} = -idx1;
                   jC                 = jC + 1;
                   
                end
                
                moveW    = blockJump(end)+1:currMD.nFrames_;
                xCorr    = circularNanCrossCor(outMap(:,moveW(1)),outMap(:,moveW(1)-1));
                [~,idx1] = max(abs(xCorr));
                
                outMap(:,moveW)    = circshift(outMap(:,moveW),-idx1);
                jump.location{jC}  = moveW;
                jump.magnitude{jC} = idx1;
                jC                 = jC + 1;
                
            else
                
                for iJ = blockJump(1):blockJump(end)%floor(numel(blockJump)/2)
                    
                    xCorr    = circularNanCrossCor(outMap(:,iJ+1),outMap(:,iJ));
                    [~,idx1] = max(xCorr);
                    
                    outMap(:,iJ+1)     = circshift(outMap(:,iJ+1),-idx1+1);
                    jump.location{jC}  = iJ;
                    jump.magnitude{jC} = -idx1+1;
                    jC                 = jC + 1;
                    
                end
                
            end
            
            
        end
    end
    
end

end