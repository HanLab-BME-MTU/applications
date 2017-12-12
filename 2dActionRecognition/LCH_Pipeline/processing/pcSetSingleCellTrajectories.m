% pcSetSingleCellTrajectories
% Input: detected cells
% Output: trajectories at the single cell level, held at the local tracking folder 

% Assaf Zaritsky, Nov 2017

function [] = pcSetSingleCellTrajectories(params,dirs)

cellTYXFname = [dirs.tracking 'cellIdTYX.mat'];

if exist(cellTYXFname,'file') && ~params.always
    fprintf(sprintf('Single cell trajectories for %s exist, finishing\n',cellTYXFname));
    return;
end

nCurCell = 0;
curCellsInds = [];
curTrajectories = [];

load([dirs.tracking 'trackingOutput.mat']); % nTrajFinal, allTrajectories, allTrajectoriesByStartTime

if ~isfield(params,'sTime')
    params.sTime = 1;
end

cellTYX = {};
for t = params.sTime : min(params.nTime - params.frameJump - 1,length(allTrajectoriesByStartTime))    
    if ~isempty(allTrajectoriesByStartTime{t})
        
        newTrajectories = allTrajectoriesByStartTime{t}.trajectories;
        nNewTrajectories = length(newTrajectories);
        
        if nNewTrajectories > 0            
            for icell = 1 : nNewTrajectories
                cellTYX{nCurCell+icell}.ts = t:(t+length(newTrajectories(icell).y)-1);
                cellTYX{nCurCell+icell}.ys = newTrajectories(icell).y;
                cellTYX{nCurCell+icell}.xs = newTrajectories(icell).x;
            end
            
            curTrajectories = [curTrajectories newTrajectories];
            curCellsInds = [curCellsInds (nCurCell+1):(nCurCell+nNewTrajectories)];
            nCurCell = nCurCell + nNewTrajectories;                                                
        end
    end
    
    for i = 1 : length(curTrajectories)
        curTraj = curTrajectories(i);
        curX = curTraj.x(t - curTraj.startFrame + 1);
        curY = curTraj.y(t - curTraj.startFrame + 1);
        
        text(curX,max(curY-50,10),sprintf(num2str(curCellsInds(i))),'color','w','FontSize',15);
    end
    
    [curTrajectories,curCellsInds] = excludeByEndFrame(curTrajectories,curCellsInds,t);
end

assert(length(cellTYX) == nCurCell);

save(cellTYXFname,'cellTYX');
end

%%
function [newTrajectories,newCellsInds] = excludeByEndFrame(trajectories,cellsInds,curFrame)
newTrajectories = [];
newCellsInds = [];
for  i = 1 : length(trajectories)
    if trajectories(i).endFrame ~= curFrame
        newTrajectories = [newTrajectories trajectories(i)];
        newCellsInds = [newCellsInds cellsInds(i)];
    end
end
end