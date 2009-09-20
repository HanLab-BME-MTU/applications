function plusTipMakeHistograms(runInfo,reclassMatrix)
% plusTipMakeHistograms saves growth, fgap, and bgap speed histograms, as well as
% a txt file containing the raw values for growth, fgap, and bgap.  
%
% reclassMatrix is the matrix where fgaps identified as continuation of growth are
% closed; that is, where trackType=5, the growth phase before, the fgap,
% and the growth phase after are merged into one subtrack, so that the
% average velocity values stored here are slightly different from that of
% projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix, where
% the rows are not merged.

homeDir=pwd;
cd(runInfo.metaDir)

histDir=[runInfo.metaDir filesep 'histograms'];
mkdir(histDir)


% here are the actual speeds according to their kinds
pop1=    reclassMatrix(reclassMatrix(:,5)==1,4);  % growth
pop2=    reclassMatrix(reclassMatrix(:,5)==2,4);  % fgap
pop3=abs(reclassMatrix(reclassMatrix(:,5)==3,4)); % bgap

% put populations into a matrix and write them into a text file
M=nan(max([length(pop1) length(pop2) length(pop3)]),3);
M(1:length(pop1),1)=pop1;
M(1:length(pop2),2)=pop2;
M(1:length(pop3),3)=pop3;
dlmwrite([runInfo.metaDir filesep 'growth_Fgap_Bgap_Distrib.txt'], M, 'precision', 3,'delimiter', '\t','newline', 'pc');

% create x-axis bins spanning all values
n=linspace(min([pop1;pop2;pop3]),max([pop1;pop2;pop3]),25);

% bin the samples
[x1,dummy] = histc(pop1,n); % growth
[x2,dummy] = histc(pop2,n); % fgap
[x3,dummy] = histc(pop3,n); % bgap

% put the binned values into a matrix for the stacked plot
M=nan(max([length(x1) length(x2) length(x3)]),3);
M(1:length(x1),1)=x1;
M(1:length(x2),2)=x2;
M(1:length(x3),3)=x3;

% make the plot
figure
bar(n,M,'stack')
colormap([1 0 0; 0 0 1; 0 1 0])
legend('growth','fgap','bgap','Location','best')
title('Stacked Speed Distributions')
xlabel('speed (um/min)');
ylabel('frequency of tracks');
saveas(gcf,[histDir filesep 'stackedHist.fig'])
saveas(gcf,[histDir filesep 'stackedHist.tif'])
close(gcf)

figure;
% growth
if ~isempty(x1)
    bar(n,x1,'r')
    title('growth speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'growthHist.fig'])
    saveas(gcf,[histDir filesep 'growthHist.tif'])
end
close(gcf)

figure
% fgap
if ~isempty(x2)
    bar(n,x2,'b')
    title('fgap speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'fgapHist.fig'])
    saveas(gcf,[histDir filesep 'fgapHist.tif'])
end
close(gcf)

figure
% bgap
if ~isempty(x3)
    bar(n,x3,'g')
    title('bgap speed distribution')
    xlabel('speed (um/min)');
    ylabel('frequency of tracks');

    saveas(gcf,[histDir filesep 'bgapHist.fig'])
    saveas(gcf,[histDir filesep 'bgapHist.tif'])
end
close(gcf)

cd(homeDir)