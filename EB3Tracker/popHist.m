function popHist(projData)

close all


% get runInfo in correct format
if nargin<1 || isempty(projData)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
    [fileName,pathName]=uigetfile('*.mat','Please select projData from META directory');
    projData=load([pathName filesep fileName]);
    projData=projData.projData;
else
    % adjust for OS
    if ~isfield(projData,'imDir') || ~isfield(projData,'anDir')
        error('--popHist: first argument should be a structure with fields imDir and anDir');
    else
        [projData.anDir] = formatPath(projData.anDir);
        [projData.imDir] = formatPath(projData.imDir);
    end
end


homeDir=pwd;
metaDir=[projData.anDir filesep 'meta'];
cd(metaDir)

% re-assign seg/gap data to shorter name
a=projData.nTrack_start_end_velMicPerMin_class_lifetime;

% get indices of segments and gaps according to their kinds
segIdx=find(a(:,5)==1);
fgapIdx=find(a(:,5)==2);
bgapIdx=find(a(:,5)==3);
ugapIdx=find(a(:,5)==4);

figure(1);
hold on

% segments
if ~isempty(a(segIdx))
    [x1 x2]=hist(a(segIdx,4),20);
    bar(x2,x1,'b')
end

% forward gaps
if ~isempty(a(fgapIdx))
    [x1 x2]=hist(a(fgapIdx,4),20);
    bar(x2,x1,'r')
end

% backward gaps
if ~isempty(a(bgapIdx))
    [x1 x2]=hist(a(bgapIdx,4),20);
    bar(x2,x1,'g')
end

% unclassified gaps
if ~isempty(a(ugapIdx))
    [x1 x2]=hist(a(ugapIdx,4),20);
    bar(x2,x1,'c')
end

xlabel('velocity (microns/minute)');
ylabel('frequency');
legend('segments (growth)','forward gaps (pause/false neg)','backward gaps (shrinkage)','unclassified','Location','best')

saveas(gcf,'compositeHist.fig')
saveas(gcf,'compositeHist.tif')

figure(2);
% segments
if ~isempty(a(segIdx))
    [x1 x2]=hist(a(segIdx,4),20);
    bar(x2,x1,'b')
    xlabel('velocity (microns/minute)');
    ylabel('frequency');
    title('histogram of average segment velocities')
    saveas(gcf,'segsHist.fig')
    frame = getframe(gca);
    [I,map] = frame2im(frame);
    imwrite(I,[pwd filesep 'segsHist.tif'],'tif')
end

figure(3);
% forward gaps
if ~isempty(a(fgapIdx))
    [x1 x2]=hist(a(fgapIdx,4),20);
    bar(x2,x1,'r')
    xlabel('velocity (microns/minute)');
    ylabel('frequency');
    title('histogram of average forward gap velocities')
    saveas(gcf,'fgapsHist.fig')
    frame = getframe(gca);
    [I,map] = frame2im(frame);
    imwrite(I,[pwd filesep 'fgapsHist.tif'],'tif')
end

figure(4);
% backward gaps
if ~isempty(a(bgapIdx))
    [x1 x2]=hist(a(bgapIdx,4),20);
    bar(x2,x1,'g')
    xlabel('velocity (microns/minute)');
    ylabel('frequency');
    title('histogram of average backward gap velocities')
    saveas(gcf,'bgapsHist.fig')
    frame = getframe(gca);
    [I,map] = frame2im(frame);
    imwrite(I,[pwd filesep 'bgapsHist.tif'],'tif')
end

figure(5);
% unclassified gaps
if ~isempty(a(ugapIdx))
    [x1 x2]=hist(a(ugapIdx,4),20);
    bar(x2,x1,'c')
    xlabel('velocity (microns/minute)');
    ylabel('frequency');
    title('histogram of average unclassified gap velocities')
    saveas(gcf,'ugapsHist.fig')
    frame = getframe(gca);
    [I,map] = frame2im(frame);
    imwrite(I,[pwd filesep 'ugapsHist.tif'],'tif')
end

figure(6);
hist(projData.pair2pairDiffPix,20)
xlabel('difference in displacement (pixels)');
ylabel('frequency');
title('difference in displacement between consecutive frame pairs')
saveas(gcf,'dispDiffHist.fig')
frame = getframe(gca);
[I,map] = frame2im(frame);
imwrite(I,[pwd filesep 'dispDiffHist.tif'],'tif')

close all
cd(homeDir)


