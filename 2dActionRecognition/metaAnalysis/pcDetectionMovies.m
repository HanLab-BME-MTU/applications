
function [] = pcDetectionMovies()

% Creates detction movies.

% Assaf Zaritsky, July 2015

warning('off','all');

always = false;

timePerFrame = 1; % minutes

FPosition = [0 0 1000 1000];
APosition = [0.1 0.1 0.9 0.9];

analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';

outdir = [analysisDirname 'Movies/detectionMovies/'];

metaDataFname = [analysisDirname 'MetaData/Experiments20151023.mat'];

load(metaDataFname);

% Create movie for each task
for i = 1 : 1 : metaData.tasks.N
    if i > metaData.tasks.N
        return;
    end
    curExp = metaData.tasks.exps(i);
    curTask = metaData.tasks.tasks(i);
    curFname = metaData.experiments.fnames{curExp};
    if curTask <= metaData.experiments.n1{curExp}
        curSource = metaData.experiments.source1{curExp};
    else
        curSource = metaData.experiments.source2{curExp};
    end
    
    %% make sure pcInitiateData was run earlier
    mdFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
        curFname '_s' sprintf('%02d',curTask) filesep...
        curFname '_s' sprintf('%02d',curTask) '.mat'];
    
    if ~exist(mdFname,'file')
        mdFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
            curFname '_s' sprintf('%d',curTask) filesep...
            curFname '_s' sprintf('%d',curTask) '.mat'];
    end
    
    if ~exist(mdFname,'file')
        continue;
    end
    
    MD =  MovieData.load(mdFname);
    
    movieFname = [outdir curFname '_s' sprintf('%d',curTask) '_detection.avi'];
    
    if exist(movieFname,'file') && ~always
        continue;
    end
    
    %     aviobj = avifile(movieFname,'fps',10,'compression','None');
    vwriter = VideoWriter(movieFname);
    vwriter.FrameRate = 10;
    open(vwriter);
    
    W = nan; H = nan;
    for j = 1 : 1000        
        detectPPDir = [analysisDirname 'Data/' curSource filesep curFname filesep...
            curFname '_s' sprintf('%02d',curTask) filesep 'detectCells/detectionsPP/'];
        
        detectPPFname = [detectPPDir sprintf('%03d',j) '_backtrack.mat'];
        
        if ~exist(detectPPFname,'file')
            continue;
        end
        
        I = MD.getChannel(1).loadImage(j);
        
        load(detectPPFname);%'detections','combinedImage'
        
        %         Ienergy = imread([detectVisDir sprintf('%03d',j) '_combinedImage.eps']);
        %         Ibacktrack = imread([detectPPVisDir sprintf('%03d',j) '_backtrack.eps']);
        
        h = figure('visible','off');%h = figure('visible','off');                
        subplot(1,2,1)
        imagesc(combinedImage); colormap(gray);
        hold on;
        plot(detections.xLowRes,detections.yLowRes,'go','MarkerSize',3,'LineWidth',1);
        plot(detections.xBacktrackLowRes,detections.yBacktrackLowRes,'co','MarkerSize',3,'LineWidth',1);
        haxes = get(h,'CurrentAxes');
        set(haxes,'XTick',[]);
        set(haxes,'XTickLabel',[]);
        set(haxes,'YTick',[]);
        set(haxes,'YTickLabel',[]);
        axis equal; axis off;        
        hold off;
        
        subplot(1,2,2)
        imagesc(I); colormap(gray);        
        hold on;
        plot(detections.x,detections.y,'go','MarkerSize',3,'LineWidth',1);
        plot(detections.xBacktrack,detections.yBacktrack,'co','MarkerSize',3,'LineWidth',1);
        text(size(I,1)-1000,size(I,2)-550,sprintf('%d minutes',round(j*timePerFrame)),'color','w','FontSize',15);
        haxes = get(h,'CurrentAxes');
        set(haxes,'XTick',[]);
        set(haxes,'XTickLabel',[]);
        set(haxes,'YTick',[]);
        set(haxes,'YTickLabel',[]);
        axis equal; axis off;
        %         axisHandle= findobj(h,'type','axes');
        %         set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off');
        
        %         axisHandle= findobj(h,'type','axes');
        %         set(axisHandle,'Position',APosition,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off');
        %         set(h,'Color','w','Position',FPosition,'PaperPositionMode','auto');
        %         set(h,'Color','w','PaperPositionMode','auto');
        hold off;
        
        drawnow; 
        pause(0.05);

        %         export_fig_biohpc('~/tmp.tif');
        %         movieFrame = imread('~/tmp.tif');
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
        %         movieFrameResized = uint8(movieFrameResized);
        movieFrame.cdata = movieFrameResized;       
        
        %         movieFrame = getframe(gcf);
        %         movieFrame = getframe(h);
        %         aviobj = addframe(aviobj, movieFrame);       
        writeVideo(vwriter,movieFrame);        
        close all;
        fprintf(sprintf('frame %d\n',j));
    end
    %     aviobj = close(aviobj);
    close(vwriter);
    
    fprintf(sprintf('done creating detection movie %d (H: %d-%d, W:%d-%d)\n,',curTask,minH,maxH,minW,maxW));
end
end