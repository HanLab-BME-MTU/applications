% update single cell DB
function [] = singleCellMovies(dbFname,analysisDname)

addpath(genpath('/home2/azaritsky/code/common'));

%% File names
if nargin == 0
    dbFname = 'singleCellLabelDB.mat';
end

if nargin < 2
    analysisDname = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/';
end

pixelSize = 0.325;

FOVRadius = 35 / pixelSize;

dbFname = [analysisDname 'SingleCellAnalysis/' dbFname];
dataDname = [analysisDname 'Data/'];

load(dbFname); % cellInfoDB,TYXDB

ncells = length(cellInfoDB);
for icell = 1 : ncells
    createSingleCellMovie(cellInfoDB{icell},analysisDname,FOVRadius);
end

end

%%
function [] = createSingleCellMovie(cellLabel,analysisDname,FOVRadius)

cellStr = [num2str(cellLabel.serialNum) '_' cellLabel.fname '_t' num2str(cellLabel.ts(1)) '_y' num2str(round(cellLabel.ys(1))) '_x' num2str(round(cellLabel.xs(1)))];

movieFname = [analysisDname filesep 'SingleCellAnalysis/Movies/' cellStr '.avi'];

if exist(movieFname,'file')
    fprintf(sprintf('%s exists\n',movieFname));
    return;
end

mdFname = [analysisDname filesep 'Data' filesep cellLabel.source filesep cellLabel.fname(1:end-4) filesep cellLabel.fname filesep cellLabel.fname '.mat']; % source, fname

MD =  MovieData.load(mdFname);

vwriter = VideoWriter(movieFname,'Uncompressed AVI');
vwriter.FrameRate = 7;
open(vwriter);

I = MD.getChannel(1).loadImage(1);
bby0 = round(max(1,cellLabel.ys(1) - FOVRadius));
bby1 = round(min(size(I,1),cellLabel.ys(1) + FOVRadius));
bbx0 = round(max(1,cellLabel.xs(1) - FOVRadius));
bbx1 = round(min(size(I,2),cellLabel.xs(1) + FOVRadius));

W = nan; H = nan;
for t = cellLabel.ts
    I = MD.getChannel(1).loadImage(t);
    
    Ibb = I(bby0:bby1,bbx0:bbx1);
    
    h = figure('visible','off'); imagesc(Ibb); colormap(gray);
    hold on;
    text(10,10,['Cell #' num2str(cellLabel.serialNum) ': ' strrep(cellLabel.fname,'_',' ') '(' num2str(cellLabel.excelID) ')'],'color','w','FontSize',15);
    text(10,20,['Labels: ' strjoin(cellLabel.labels,', ')],'color','w','FontSize',12);
    text(10,30,['Descriptor: ' cellLabel.desc],'color','w','FontSize',12);
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