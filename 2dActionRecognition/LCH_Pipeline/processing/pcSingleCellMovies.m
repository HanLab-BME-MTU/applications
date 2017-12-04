
function [] = pcSingleCellMovies(MD,params,dirs)

% Create movie for every single cell at /project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Cells/Movies

% Assaf Zaritsky, June 2016

outdir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Cells/Movies/';

cellTYXFname = [dirs.tracking 'cellIdTYX.mat'];
load(cellTYXFname);% cellTYX

nCells = length(cellTYX);

FOVRadius = 35 / params.pixelSize;

for icell = 1 : nCells
    ts = cellTYX{icell}.ts;
    xs = cellTYX{icell}.xs;
    ys = cellTYX{icell}.ys;
    createSingleCellMovie(MD,xs,ys,ts,FOVRadius,outdir,dirs.expname,icell);
end
end

%%
function [] = createSingleCellMovie(MD,xs,ys,ts,FOVRadius,outdir,expname,icell)
movieDname = [outdir filesep expname filesep];

if ~exist(movieDname,'dir')
    unix(sprintf('mkdir %s',movieDname));
end

movieFname = [movieDname filesep num2str(icell) '.avi']; % '_t' num2str(cellLabel.ts(1)) '_y' num2str(round(cellLabel.ys(1))) '_x' num2str(round(cellLabel.xs(1)))

if exist(movieFname,'file')
    fprintf(sprintf('%s exists\n',movieFname));
    return;
end

vwriter = VideoWriter(movieFname,'Uncompressed AVI');
vwriter.FrameRate = 7;
open(vwriter);

I = MD.getChannel(1).loadImage(1);
bby0 = round(max(1,ys(1) - FOVRadius));
bby1 = round(min(size(I,1),ys(1) + FOVRadius));
bbx0 = round(max(1,xs(1) - FOVRadius));
bbx1 = round(min(size(I,2),xs(1) + FOVRadius));

W = nan; H = nan;
for t = ts
    I = MD.getChannel(1).loadImage(t);
    
    Ibb = I(bby0:bby1,bbx0:bbx1);
    
    h = figure('visible','off'); imagesc(Ibb); colormap(gray);
    hold on;
    %     text(10,10,['Cell #' num2str(cellLabel.serialNum) ': ' strrep(cellLabel.fname,'_',' ') '(' num2str(cellLabel.excelID) ')'],'color','w','FontSize',15);
    %     text(10,20,['Labels: ' strjoin(cellLabel.labels,', ')],'color','w','FontSize',12);
    %     text(10,30,['Descriptor: ' cellLabel.desc],'color','w','FontSize',12);
    timePerFrame = 1;
    text(bbx1-bbx0-40,bby1-bby0-10,sprintf('%d minutes',round(t*timePerFrame)),'color','w','FontSize',12);
    
    haxes = get(h,'CurrentAxes');
    set(haxes,'XTick',[]);
    set(haxes,'XTickLabel',[]);
    set(haxes,'YTick',[]);
    set(haxes,'YTickLabel',[]);
    
    hold off;
    
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
end

end