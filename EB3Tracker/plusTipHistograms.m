function plusTipHistograms(projData)
% plusTipHistograms saves growth, pause, and shrinkage speed histograms and an
% Excel spreadsheet containing those three speed distribution




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
metaDir=[formatPath(projData.anDir) filesep 'meta'];
cd(metaDir)

histDir=[metaDir filesep 'histograms'];
if ~isdir(histDir)
    mkdir(histDir)
end

% re-assign seg/gap data to shorter name
a=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;

% get indices of segments and gaps according to their kinds
segIdx=find(a(:,5)==1);
fgapIdx=find(a(:,5)==2);
bgapIdx=find(a(:,5)==3);

% here are the actual speeds
pop1=a(segIdx,4);
pop2=a(fgapIdx,4);
pop3=abs(a(bgapIdx,4));

% write data into Excel spreadsheet
% if ispc
%     letters={'A' 'B' 'C' 'D'};
%     names{1,1}='growth';
%     names{2,1}='pause';
%     names{3,1}='shrinkage';
%     for iName=1:length(names)
%         switch iName
%             case 1
%                 tempMat=pop1;
%                 % put headers at the top of the file
%                 range=[letters{iName+1} num2str(1) ':' letters{iName+3} num2str(1)];
%                 xlswrite([metaDir filesep 'speedDistributions'],{'growth', 'pause' 'shrinkage'},range);
%             case 2
%                 tempMat=pop2;
%             case 3
%                 tempMat=pop3;
%         end
%         % write the growth, pause, and shrinkage distributions into the file
%         [r c]=size(tempMat);
%         range=[letters{iName+1} num2str(3) ':' letters{iName+1} num2str(3+r-1)];
%         xlswrite([metaDir filesep 'speedDistributions'],tempMat,range);
%     end
% else
M=nan(max([length(pop1) length(pop2) length(pop3)]),3);
M(1:length(pop1),1)=pop1;
M(1:length(pop2),2)=pop2;
M(1:length(pop3),3)=pop3;
dlmwrite([metaDir filesep 'growthPauseShrinkSpeedDistributions.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');
% end

% create x-axis bins spanning all costs in sample
n=linspace(min([pop1;pop2;pop3]),max([pop1;pop2;pop3]),25);

% bin the samples
[x1,nbins1] = histc(pop1,n); % forward
[x2,nbins2] = histc(pop2,n); % backward
[x3,nbins3] = histc(pop3,n); % backward

M=nan(max([length(x1) length(x2) length(x3)]),3);
M(1:length(x1),1)=x1;
M(1:length(x2),2)=x2;
M(1:length(x3),3)=x3;

% make the plot
figure
bar(n,M,'stack')
colormap([1 0 0; 0 0 1; 0 1 0])
legend('growth','pause','shrinkage','Location','best')
title('Sub-Track Speed Distribution')
xlabel('speed (um/min)');
ylabel('frequency of tracks');
saveas(gcf,[histDir filesep 'stackedHist.fig'])
saveas(gcf,[histDir filesep 'stackedHist.tif'])

close(gcf)

figure;
% growth
if ~isempty(pop1)
    [x1 x2]=hist(pop1,25);
    bar(x2,x1,'b')
    title('growth speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'growthHist.fig'])
    saveas(gcf,[histDir filesep 'growthHist.tif'])
end
close(gcf)

figure
% pause
if ~isempty(pop2)
    [x1 x2]=hist(pop2,25);
    bar(x2,x1,'r')
    title('pause speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'pauseHist.fig'])
    saveas(gcf,[histDir filesep 'pauseHist.tif'])
end
close(gcf)

figure
% shrinkage
if ~isempty(pop3)
    [x1 x2]=hist(pop3,25);
    bar(x2,x1,'g')
    title('shrinkage speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'shrinkageHist.fig'])
    saveas(gcf,[histDir filesep 'shrinkageHist.tif'])
end
close(gcf)

cd(homeDir)


