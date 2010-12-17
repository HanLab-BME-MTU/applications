function [groupCorr]=groupCorrResults(groupCorr,corr)
[grpBin,grpFrame]=size(groupCorr);
[sglBin,sglFrame]=size(     corr);

maxBin  =max(grpBin  ,sglBin);
maxFrame=max(grpFrame,sglFrame);

for bin=1:maxBin
    for frame=1:maxFrame
        if bin>sglBin || frame>sglFrame || isempty(corr(bin,frame))
            corr=initCorr(corr,bin,frame);
        end
        if bin>grpBin || frame>grpFrame || isempty(groupCorr(bin,frame))
            groupCorr=initCorr(groupCorr,bin,frame);
        end        
        
        grpNumDr=length(groupCorr(bin,frame).RMean);
        newNumDr=length(     corr(bin,frame).RMean);
        
        groupCorr(bin,frame).RMean          = vertcat(padarray(groupCorr(bin,frame).RMean   , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).RMean    , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).RSTD           = vertcat(padarray(groupCorr(bin,frame).RSTD    , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).RSTD     , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).funcVel        = vertcat(padarray(groupCorr(bin,frame).funcVel , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).funcVel  , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).funcVelx       = vertcat(padarray(groupCorr(bin,frame).funcVelx, [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).funcVelx , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).funcVely       = vertcat(padarray(groupCorr(bin,frame).funcVely, [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).funcVely , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).Ntot           = vertcat(padarray(groupCorr(bin,frame).Ntot    , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).Ntot     , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        % the angle correlation:
        groupCorr(bin,frame).RMean4CosaMean = vertcat(padarray(groupCorr(bin,frame).RMean4CosaMean, [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).RMean4CosaMean, [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).RMean4CosaSTD  = vertcat(padarray(groupCorr(bin,frame).RMean4CosaSTD , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).RMean4CosaSTD , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).cosaMean       = vertcat(padarray(groupCorr(bin,frame).cosaMean      , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).cosaMean      , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        groupCorr(bin,frame).cosaSTD        = vertcat(padarray(groupCorr(bin,frame).cosaSTD       , [0,max(grpNumDr,newNumDr)-grpNumDr], NaN  ,'post'), padarray(     corr(bin,frame).cosaSTD       , [0,max(grpNumDr,newNumDr)-newNumDr], NaN  ,'post'));
        % the in-sheet positions:
        groupCorr(bin,frame).c2edMean       = vertcat(groupCorr(bin,frame).c2edMean,corr(bin,frame).c2edMean);
        groupCorr(bin,frame).c2edSTD        = vertcat(groupCorr(bin,frame).c2edSTD ,corr(bin,frame).c2edSTD);
        
        
        
        % also store the full vectors:
        grpR=groupCorr(bin,frame).R;
        sglR=corr(bin,frame).R;
        
        grpc2cBin=length(grpR);
        sglc2cBin=length(sglR);        
        
        maxc2cBin=max(grpc2cBin,sglc2cBin);
        for c2cBin=1:maxc2cBin
            if c2cBin>sglc2cBin || isempty(sglR(c2cBin))
                sglR=initR(sglR,c2cBin);
            end
            if c2cBin>grpc2cBin || isempty(grpR(c2cBin))
                grpR=initR(grpR,c2cBin);
            end           
            
            grpR(c2cBin).pts2ptR = vertcat(grpR(c2cBin).pts2ptR,sglR(c2cBin).pts2ptR);
            grpR(c2cBin).xVelPts = vertcat(grpR(c2cBin).xVelPts,sglR(c2cBin).xVelPts);
            grpR(c2cBin).yVelPts = vertcat(grpR(c2cBin).yVelPts,sglR(c2cBin).yVelPts);
            grpR(c2cBin).xVelCtr = vertcat(grpR(c2cBin).xVelCtr,sglR(c2cBin).xVelCtr);
            grpR(c2cBin).yVelCtr = vertcat(grpR(c2cBin).yVelCtr,sglR(c2cBin).yVelCtr);
            grpR(c2cBin).cosda   = vertcat(grpR(c2cBin).cosda  ,sglR(c2cBin).cosda);
            % increase the total number of correlation pairs:
            grpR(c2cBin).Ntot    = vertcat(grpR(c2cBin).Ntot,sglR(c2cBin).Ntot);
        end
        groupCorr(bin,frame).R=grpR;
    end
end

function newCorr=initCorr(newCorr,newBin,newFrame)
newCorr(newBin,newFrame).RMean          = [];
newCorr(newBin,newFrame).RSTD           = [];
newCorr(newBin,newFrame).funcVel        = [];
newCorr(newBin,newFrame).funcVelx       = [];
newCorr(newBin,newFrame).funcVely       = [];
newCorr(newBin,newFrame).Ntot           = [];
% the angle correlation:
newCorr(newBin,newFrame).RMean4CosaMean = [];
newCorr(newBin,newFrame).RMean4CosaSTD  = [];
newCorr(newBin,newFrame).cosaMean       = [];
newCorr(newBin,newFrame).cosaSTD        = [];
% the in-sheet positions:
newCorr(newBin,newFrame).c2edMean       = [];
newCorr(newBin,newFrame).c2edSTD        = [];

% the full vectors:
newCorr(newBin,newFrame).R=[];


function newR=initR(newR,newRBin)
newR(newRBin).pts2ptR = [];
newR(newRBin).xVelPts = [];
newR(newRBin).yVelPts = [];
newR(newRBin).xVelCtr = [];
newR(newRBin).yVelCtr = [];
newR(newRBin).cosda   = [];
% increase the total number of correlation pairs:
newR(newRBin).Ntot    = [];