%Script for subplot all superresolution data
% Script to plot the p-value figures as an array: rows: dissociation rate,
% columns: different receptor density. It is divided into two parts: first the
% dynamic data and then the super-resolution.

saveDir='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170525/superRes/results/figures';

rDtarget = {'rD4'};%,'rD8','rD10','rD12','rD14','rD16'};
aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
dRT={'dR0p5'};
% values of rD, aP and lR of probe
rDvals = 4;
aPvals = 0.2:0.1:0.8;
dRVals=0.25:0.25:2;
dRdir = {'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};
lRtargetSR ={'lR1'};%,'lR0p3','lR0p4','lR0p5'};
lRDir={'lR0p91';'lR0p93';'lR0p95';'lR0p97';'lR0p99';'lR1'};

%%%%%%%%%%%%%%%%%%%%%%%%%%super-resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsDirectorySR='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170525/superRes/results';

%
dRindx=1;

%%%%%%%%%%%%% plot list
%1- plotIndx=[1,7,13,19,25,31,37,43]
%2- plotIndx=[2,8,14,20,26,32,38,44]
%3- plotIndx=[3,9,15,21,27,33,39,45]
%4- plotIndx=[4,10,16,22,28,34,40,46]
%5- plotIndx=[5,11,17,23,29,35,41,47]
%6- plotIndx=[6,12,18,24,30,36,42,48]




    %     for dRindx=1:length(dRdir)
    temp= load([resultsDirectorySR,filesep,rDtarget{1},dRT{1},aPtarget{1},lRtargetSR{1 },filesep,'pMatrix_',lRtargetSR{1},'.mat']);
    for lRIndex=1:length(lRDir)
    pMatrix=temp.pMatrix;
    pMatrix(pMatrix<0.05)=0.05;
    subplot(5,6,lRIndex);
    ax = gca;
    load('MyColormaps','mycmap')
    colormap(ax,mycmap)
    imagesc(aPvals,dRVals,pMatrix(:,:,1,lRIndex));
%     %to create the grid
     line(repmat([0 2],9,1).',repmat(0.125:0.25:2.125,2,1),'Color',ones(1,3)*0.2)
     line(repmat(0.15:0.1:0.85,2,1),repmat([0; 2.5],1,8),'Color',ones(1,3)*0.2)
    
    caxis([0.05 1])
    % configurations to have the graphic ploted in the "normal" direction.
    axH = gca; %ax = gca returns the handle to the current axes for the
    % current figure. If an axes does not exist, then gca creates an axes and
    % returns its handle. You can use the axes handle to query and modify axes
    % properties.
    set(axH,'YDir','normal');
    
    %%%%%%%%%%%%%%%removing number in the axis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%                    set(gca,'YTick',[])
%                     set(gca,'XTick',[])
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dRindx=dRindx+1;
    end 

% end
% %include color bar
   ax=gca;
   pos=get(gca,'pos');
  set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
   pos=get(gca,'pos');
   hc=colorbar('location','northoutside','position',[pos(1) pos(2)+pos(4) pos(3) 0.03]);
   set(hc,'xaxisloc','top');
% %
% %save figure
%
% %Save figure
figH = gcf;
set(figH,'Name',lRtargetSR{1});
resDir=[saveDir,filesep,filesep,rDtarget{1},aPtarget{1},lRtargetSR{1 }];
mkdir(resDir)
outFile =[resDir ,filesep,'pMatrixPlot']; % %
%                    saveas(figH,outFile,'png');
  saveas(figH,outFile,'fig');
% %
%
%
%
% %inclues a title for

% [ax,h]=subtitle([rDtarget{1},aPtarget{1},lRtarget{1},dRT{1}]);
