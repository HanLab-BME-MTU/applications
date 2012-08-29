function analyzekMTs(MD)

close all
gfpIndx = find(strcmpi('gfp',{MD.channels_.fluorophore_}));
mcherryIndx = find(strcmpi('mcherry',{MD.channels_.fluorophore_}));
% Get the threshold indices
iThresProc = MD.getProcessIndex('ThresholdProcess',1,0);
iMaskRefProc = MD.getProcessIndex('MaskRefinementProcess',1,0);
iBackMaskProc = MD.getProcessIndex('BackgroundMasksProcess',1,0);
iKinDetProc = MD.getProcessIndex('PointSourceDetectionProcess',1,0);
iTrackProc = MD.getProcessIndex('TrackingProcess',Inf,0);
iKinTrackProc = iTrackProc(arrayfun(@(x)isequal(MD.processes_{x}.funParams_.DetProcessIndex,iKinDetProc),...
    iTrackProc));
iCometDetProc = MD.getProcessIndex('CometDetectionProcess',1,0);

% Filter tracks by lifetime
tracks = MD.processes_{iKinTrackProc}.loadChannelOutput(1);
tracksLifetime=arrayfun(@(x) diff(x.seqOfEvents(:,1)),tracks);
kinTracks=tracks(tracksLifetime>.8*MD.nFrames_);
nTracks=numel(kinTracks);

% Get comet info
movieInfo = MD.processes_{iCometDetProc}.loadChannelOutput(2);
allCometPos = arrayfun(@(x) [x.xCoord(:,1) x.yCoord(:,1)],movieInfo,'Unif',0);

%% Calculat closest comet
tic
xKin=cell(nTracks,1);
yKin=cell(nTracks,1);
xkMT=cell(nTracks,1);
ykMT=cell(nTracks,1);
AkMT=cell(nTracks,1);
IkMT=cell(nTracks,1);
threshold=8;

for iTrack=1:nTracks;
    xKin{iTrack}=kinTracks(iTrack).tracksCoordAmpCG(1:8:end);
    yKin{iTrack}=kinTracks(iTrack).tracksCoordAmpCG(2:8:end);
    
    nTrackedFrames =numel(xKin{iTrack});
    xkMT{iTrack}=NaN(nTrackedFrames,1);
    ykMT{iTrack}=NaN(nTrackedFrames,1);
    AkMT{iTrack}=NaN(nTrackedFrames,1);
    IkMT{iTrack}=NaN(nTrackedFrames,1);
    for t=1:nTrackedFrames
        idx=KDTreeBallQuery(allCometPos{t},[xKin{iTrack}(t) yKin{iTrack}(t)],threshold);
        if ~isempty(idx{1})
            xkMT{iTrack}(t)=allCometPos{t}(idx{1}(1),1);
            ykMT{iTrack}(t)=allCometPos{t}(idx{1}(1),2);
            AkMT{iTrack}(t)=movieInfo(t).amp(idx{1}(1));
            IkMT{iTrack}(t)=movieInfo(t).maxInt(idx{1}(1));
        end
    end
end
toc


%% Create output directory

outputPath = [MD.outputDirectory_ filesep 'kMTs'];
mkClrDir(outputPath);

%% Figure with all tracks

% figure('Visible','off'); 
% MD.channels_.draw(1);
% hold on;
% for iTrack=1:nTracks;
%     plot(xKin{iTrack},yKin{iTrack},'-k');
%     plot3(xkMT{iTrack},ykMT{iTrack},AkMT{iTrack},'-r');
% end
% print('-dpng', '-r300', [outputPath filesep 'Tracks.png']); % the -r300 option is only needed when figure has non-vectorial content


% %% Movie
% 
% y0=113;
% x0=156;
% % y0=167;
% % x0=193;
% length=40;
% height=15;
% cropROI=[x0 y0 length-1 height-1];
% hKin=-1*ones(nTracks,1);
% hKinTrack=-1;
% hkMT=-1*ones(2*nTracks,1);
% dragTail=10;
% nColors=64;
% colors=jet(nColors);
% Amax = cellfun(@max,AkMT);
% ext = '.png';
% 
% zoom = 1; % zoom level compared to screen resolution
% % Configure figure
% nx = length*20;
% ny=height*20;
% h = figure('Visible', 'on', 'Position', [500 500 nx ny]);
% iptsetpref('ImshowBorder','tight');
% set(h, 'InvertHardcopy', 'off');
% set(h, 'PaperUnits', 'Points');
% set(h, 'PaperSize', [nx ny]);
% set(h, 'PaperPosition', [0 0 nx ny]); % very important
% axes('Position',[0 0 1 1]);
% axis([x0 x0+length-1 y0 y0+height-1]);
% % set(h,'DefaultLineLineSmoothing','on');
% % set(h,'DefaultPatchLineSmoothing','on');
% 
% framesPath = [MD.outputDirectory_ filesep 'frames' filesep];
% mkdir(framesPath);
% fmt = ['%0' num2str(ceil(log10(MD.nFrames_))) 'd'];
% 
% frames=10:40;
% for i=1:numel(frames)
%     t=frames(i);
%     delete(hKin(ishandle(hKin)));
%     delete(hKinTrack(ishandle(hKinTrack)));
%     delete(hkMT(ishandle(hkMT)));
%     MD.channels_.draw(t,'hAxes',gca);
%     set(gca,'XLim',[x0 x0+length-1],'YLim',[y0 y0+height-1]);
%     hold on;
%     for iTrack=1:nTracks
%         hKinTrack(iTrack,1)=plot(xKin{iTrack}(max(1,t-20):t),yKin{iTrack}(max(1,t-20):t),'-r');
%         hKin(iTrack,1)=plot(xKin{iTrack}(t),yKin{iTrack}(t),'ok','MarkerFaceColor','r');
%         colorIndx = round(AkMT{iTrack}(t)/Amax(iTrack)*64);
%         if ~isnan(colorIndx)
%             hkMT(nTracks+iTrack,1)=plot(xkMT{iTrack}(t),ykMT{iTrack}(t),'o',...
%                 'MarkerSize',4+round(AkMT{iTrack}(t)/Amax(iTrack)*20));
% %                 'MarkerSize',colors(colorIndx,:));
%         end
%     end
%     print(h, '-dpng', '-loose', ['-r' num2str(zoom*72)], [framesPath 'frame' num2str(i, fmt) ext]);
% 
%     
% end
% 
% 
% fr = '8';
% cmd = ['ffmpeg -y -r ' fr ' -i ' regexprep(framesPath,' ','\\ ') 'frame' fmt ext...
%     ' -r ' fr ' -b 50000k -bt 20000k ' regexprep(framesPath,' ','\\ ') 'movie.mov > /dev/null'];
% 
% system(cmd);
%  
%% Plot position/velocity/amplitude
% define small and large fonts
tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
sfont = {'FontName', 'Helvetica', 'FontSize', 18};
lfont = {'FontName', 'Helvetica', 'FontSize', 22};

MD.timeInterval_=1;
for iTrack=1:nTracks
    
    dx=diff(sgolayfilt(xKin{iTrack},1,11));    
    amp=IkMT{iTrack}(1:end-1);
    amp(amp==0)=NaN;
    time=(1:numel(amp))*MD.timeInterval_;

    f=figure('Visible','off','Position',[0 0 800 600],'PaperPositionMode', 'auto'); % enable resizing
    [ax h(1) h(2)]=plotyy(time,dx*MD.pixelSize_/1000/MD.timeInterval_*60,...
        time,amp/max(amp(:)));
    xlabel('Time (s)',lfont{:});
    set(h(1),'Linewidth',2,'Color','r');
    set(h(2),'Linewidth',2,'Color','b');
    set(ax(1),'YColor','r');
    set(ax(2),'YColor','b');
    set(ax, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    set(get(ax(1),'Ylabel'),'String','Kinetochore speed (microns/min)',lfont{:});
    set(get(ax(2),'Ylabel'),'String','Comet maximal intensity',lfont{:});

    print('-dpng', '-r300', [outputPath filesep 'Track' num2str(iTrack) ...
        '-Speed-intensity.png']); % the -r300 option is only needed when figure has non-vectorial content
%     print('-depsc2', '-r300', 'Speed-intensity.png'); % the -r300 option is only needed when figure has non-vectorial content
    % default density is 72, i.e., rasterize at screen resolution
%     !convert -quiet -colorspace rgb -alpha off -depth 8 -density 300 -compress LZW Speed-intensity.eps Speed-intensity.tif
    close(f)
 
    f=figure('Visible','off','Position',[0 0 800 600],'PaperPositionMode', 'auto'); % enable resizing
    [ax h(1) h(2)]=plotyy(time,xKin{iTrack}(1:end-1)*MD.pixelSize_/1000,...
        time,amp/max(amp(:)));
    xlabel('Time (s)',lfont{:});
    set(h(1),'Linewidth',2,'Color','r');
    set(h(2),'Linewidth',2,'Color','b');
    set(ax(1),'YColor','r');
    set(ax(2),'YColor','b');
    set(ax, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    set(get(ax(1),'Ylabel'),'String','Kinetochore position (microns)',lfont{:});
    set(get(ax(2),'Ylabel'),'String','Comet maximal intensity',lfont{:});
    
    print('-dpng', '-r300', [outputPath filesep 'Track' num2str(iTrack) ...
        '-Position-intensity.png']); % the -r300 option is only needed when figure has non-vectorial content
    close(f)
    % default density is 72, i.e., rasterize at screen resolution
    %     !convert -quiet -colorspace rgb -alpha off -depth 8 -density 300 -compress LZW Position-intensity.eps Position-intensity.tif
    
    %
    % dx=cellfun(@diff,xKin,'Unif',0);
    % dx=vertcat(dx{:});
    %
    % dy=cellfun(@diff,yKin,'Unif',0);
    % dy=vertcat(dy{:});
    %
    % v=(dx.^2+dy.^2).^(.5);
    % amp=horzcat(AkMT{:});
    % amp=amp(1:end-1,:)';
    % figure;
    % plot(v,amp,'.r');
    
    
    %% Plot cross correlation
%     X=xKin{iTrack}(1:end-1);
    X=dx;
    Y=amp';
    X(detectOutliers(X,5)) = NaN;
    Y(detectOutliers(Y,5)) = NaN;
    
    % Percentage of NaN
    [X,range1] = removeMeanTrendNaN(X');
    [Y,range2] = removeMeanTrendNaN(Y');
    
    [~,range1,range2] = intersect(range1{1},range2{1});
    ccL    = numel(range1);
    nLags = round(ccL/4);
    [corrFun,lags,bounds] =...
        crosscorr(X{1}(range1),Y{1}(range2),nLags);
    
    
    [~,indx]=max(corrFun);
    imin=max(1,indx-2);
    imax=min(indx+2,numel(lags));
    [coeffs]=polyfit(lags(imin:imax),corrFun(imin:imax),2);
    x0=-coeffs(2)/(2*coeffs(1));
    
    f=figure('Visible','off','Position',[0 0 800 600],'PaperPositionMode', 'auto'); % enable resizing
    plot(lags,corrFun,'Linewidth',2);
    hold on;
    plot(lags,repmat(bounds(1),1,numel(lags)),'r','Linewidth',2);
    plot(lags,repmat(bounds(2),1,numel(lags)),'g','Linewidth',2);
    ylim=get(gca,'YLim');
    plot([x0 x0],ylim,'--k','Linewidth',2);
    xlabel('Time (s)',lfont{:})
    ylabel('Correlation function',lfont{:});
    set(gca, 'LineWidth', 1.5, sfont{:}, 'Layer', 'top');
    
    
    print('-dpng', '-r300', [outputPath filesep 'Track' num2str(iTrack) ...
        'crosscorrelation.png']); % the -r300 option is only needed when figure has non-vectorial content
    % default density is 72, i.e., rasterize at screen resolution
    %     !convert -quiet -colorspace rgb -alpha off -depth 8 -density 300 -compress LZW crosscorrelation.eps crosscorrelation.tif
    close (f)
end