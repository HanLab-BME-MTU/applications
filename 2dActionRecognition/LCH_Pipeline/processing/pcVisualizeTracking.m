function [] = pcVisualizeTracking(MD,params,dirs)

close all;

trackingVisFname = [dirs.trackingVis dirs.expname '_tracking.avi'];

load([dirs.tracking 'trackingOutput.mat']); % nTrajFinal, allTrajectories, allTrajectoriesByStartTime

fprintf('visualize tracking\n');  

aviobj = avifile(trackingVisFname,'fps',10,'compression','None');

nCurCell = 0;
curCellsInds = [];
curTrajectories = [];
for t = 1 : params.nTime - params.frameJump - 1
    I = MD.getChannel(1).loadImage(t);
    
    if ~isempty(allTrajectoriesByStartTime{t})
        
        newTrajectories = allTrajectoriesByStartTime{t}.trajectories;
        nNewTrajectories = length(newTrajectories);
        
        if nNewTrajectories > 0
            curTrajectories = [curTrajectories newTrajectories];
            curCellsInds = [curCellsInds (nCurCell+1):(nCurCell+nNewTrajectories)];
            nCurCell = nCurCell + nNewTrajectories;
        end
    end
    
    
    h = figure('visible','off'); imagesc(I); colormap(gray);
    %         h = figure(); imagesc(I); colormap(gray);
    text(size(I,1)-300,size(I,2)-500,sprintf('%d minutes',round(t*params.timePerFrame)),'color','w','FontSize',15);
    for i = 1 : length(curTrajectories)
        curTraj = curTrajectories(i);
        curX = curTraj.x(t - curTraj.startFrame + 1);
        curY = curTraj.y(t - curTraj.startFrame + 1);
        text(curX,max(curY-50,10),sprintf(num2str(curCellsInds(i))),'color','w','FontSize',15);
    end
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    drawnow; pause(0.01);
    movieFrame = getframe(h);
    aviobj = addframe(aviobj, movieFrame);
    close all;
    
    [curTrajectories,curCellsInds] = excludeByEndFrame(curTrajectories,curCellsInds,t);
    
end
aviobj = close(aviobj);
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