
function [] = pcTrackingMovie(MD,params,dirs)

% Save tracking movies to /project/cellbiology/gdanuser/melanomaModel/Analysis/Movies/trackingMovies
% Note: over 2 minutes per frame from getFrame!!

% Assaf Zaritsky, November 2015 
% Nov. 2017 - redundantly redifining the cells and all here the same way as
% within pcSetSingleCellTrajectories. This is because the trajectories are
% orginized according to their time frame which is very convinient for
% visualization. But I did want the modularity to have one function with
% no I/O that does the trajectories definition. 

outdir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Movies/trackingMovies';
% '/project/cellbiology/gdanuser/melanomaModel/Analysis/Movies/trackingMovies/';

movieFname = [outdir filesep dirs.expname '_tracking.avi'];
cellTYXFname = [dirs.tracking 'cellIdTYX.mat'];

if exist(movieFname,'file') && exist(cellTYXFname,'file') && ~params.always
    fprintf(sprintf('Tracking movie %s exists, finishing\n',movieFname));
    return;
end

vwriter = VideoWriter(movieFname,'Uncompressed AVI');
vwriter.FrameRate = 10;
open(vwriter);

nCurCell = 0;
curCellsInds = [];
curTrajectories = [];

load([dirs.tracking 'trackingOutput.mat']); % nTrajFinal, allTrajectories, allTrajectoriesByStartTime

if ~isfield(params,'sTime')
    params.sTime = 1;
end

cellTYX = {};
W = nan; H = nan;
for t = params.sTime : min(params.nTime - params.frameJump - 1,length(allTrajectoriesByStartTime))
    I = MD.getChannel(1).loadImage(t);
    
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
        
    if isnan(W)
        [H,W,~] = size(movieFrame.cdata);
        minH = H;
        maxH = H;
        minW = W;
        maxW = W;
    end
    
    if H ~= size(movieFrame.cdata,1) || W ~= size(movieFrame.cdata,2)
        minH = min(H,size(movieFrame.cdata,1));
        maxH = max(H,size(movieFrame.cdata,2));
        minW = min(W,size(movieFrame.cdata,1));
        maxW = max(W,size(movieFrame.cdata,2));
    end
    
    movieFrameResized = uint8(zeros(H,W,3));
    movieFrameResized(:,:,1) = imresize(movieFrame.cdata(:,:,1),[H,W]);
    movieFrameResized(:,:,2) = imresize(movieFrame.cdata(:,:,2),[H,W]);
    movieFrameResized(:,:,3) = imresize(movieFrame.cdata(:,:,3),[H,W]);
    movieFrame.cdata = movieFrameResized;
    
    writeVideo(vwriter,movieFrame);
    close all;
    fprintf(sprintf('frame %d\n',t));
    
    close all;
    
    [curTrajectories,curCellsInds] = excludeByEndFrame(curTrajectories,curCellsInds,t);
    
end

assert(length(cellTYX) == nCurCell);

% save(cellTYXFname,'cellTYX');
close(vwriter);

fprintf(sprintf('Done creating tracking movie (H: %d-%d, W:%d-%d)\n,',minH,maxH,minW,maxW));
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