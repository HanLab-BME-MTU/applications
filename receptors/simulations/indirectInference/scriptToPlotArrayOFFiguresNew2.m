% Script to plot the p-value figures as an array: rows: dissociation rate,
% columns: labeled fraction. It is divided into two parts: first the
% dynamic data and then the super-resolution.

saveDir= '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170508/results/figures/20170518';
%%%%%%%%%%%%%%%%%%Dynamic data path%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultsDirectory='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170508/results/dyanamicData';
%%%%%%%%%%%%%%%%%%static data path%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultsDirectorySR='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170508/results/staticData';


% % name of the target
rDtarget = {'rD4'};%,'rD8','rD10','rD12','rD14','rD16'};
aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
lRtarget ={'lR0p2'};%,'lR0p3','lR0p4','lR0p5'};
dRT={'dR0p5'};
lRT={'lR1'};
% values of rD, aP and lR of probe
dRvals =0.25:0.25:2;
aPvals = 0.2:0.1:0.8;
lRStr = {'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5','lR0p6'};%,,'lR0p93','lR0p95','lR0p97','lR0p99'
%  dRdir ={'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};%{'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};
rDDir={'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'};
% calculate length of lR because it will be used for calculation of
% position in the rows
lengthLr=length(lRStr);
lengthDr=length(rDDir);
%figures

% first we need to vary in the dissociation rate, because it will be the
% variations in the number of rows
rDindx=1;

% for lRTIndx = 1 : length(lRtarget)
    
    %Iterate through association probability values per density
    for aPTIndx = 1 : length(aPtarget)
        %
        %
        %     %iterate through the different labeling ratios
%         for rDTIndx = 1 : length(rDDir)
            %
            %
            %
            %
            %
            for rowIndex=[50,50-lengthLr-1,50-2*lengthLr-2,50-3*lengthLr-3,50-4*lengthLr-4,50-5*lengthLr-5,50-6*lengthLr-6,50-7*lengthLr-7]
                temp= load([resultsDirectory,filesep,rDtarget{1},dRT{1},aPtarget{1},lRtarget{1},filesep,'pMatrix.mat']);
                pMatrix=temp.pMatrix;
                pMatrix(pMatrix<0.05)=0.05;
                %
                %
                %
                for lRindx=1:lengthLr-1
                    %                 subplot(8,7,lengthLr*dRindx)
                    subplot(8,7,rowIndex+lRindx)
                    %
                    %
                    ax = gca;
                    load('MyColormaps','mycmap')
                    colormap(ax,mycmap)
                    imagesc(aPvals,dRvals,pMatrix(:,:,rDindx,lRindx));
                    %                 %to create the grid
                    line(repmat([0 2],9,1).',repmat(0.125:0.25:2.125,2,1),'Color',ones(1,3)*0.2)
                    line(repmat(0.15:0.1:0.85,2,1),repmat([0; 2.5],1,8),'Color',ones(1,3)*0.2)
                    %
                    caxis([0.05 1])
                    %                 % configurations to have the graphic ploted in the "normal" direction.
                    axH = gca; %ax = gca returns the handle to the current axes for the
                    %                 % current figure. If an axes does not exist, then gca creates an axes and
                    %                 % returns its handle. You can use the axes handle to query and modify axes
                    %                 % properties.
                    set(axH,'YDir','normal');
                    %
                    %                 %%%%%%%%%%%%%%%removing number in the axis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % %                                  set(gca,'YTick',[])
                    % %                                  set(gca,'XTick',[])
                    %
                    %
                    %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    %
                end
                %             %increase dRindx
                rDindx=rDindx+1;
            end
            %         %increase
%         end
        %     %               colorbar
    end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%super-resolution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%
 rDindxN=8;
lRtargetSR ={'lR1'};%,'lR0p3','lR0p4','lR0p5'};
% %

for rDindx=1:length(rDDir)
    temp= load([resultsDirectorySR,filesep,rDtarget{1},dRT{1},aPtarget{1},lRtargetSR{1},filesep,'pMatrix_',lRtargetSR{1},'.mat']);
    pMatrix=temp.pMatrix;
    pMatrix(pMatrix<0.05)=0.05;
    subplot(8,7,7*rDindxN);
    ax = gca;
    load('MyColormaps','mycmap')
    colormap(ax,mycmap)
    imagesc(aPvals,dRvals,pMatrix(:,:,rDindx));
    % %to create the grid
    line(repmat([0 2],9,1).',repmat(0.125:0.25:2.125,2,1),'Color',ones(1,3)*0.2)
    line(repmat(0.15:0.1:0.85,2,1),repmat([0; 2.5],1,8),'Color',ones(1,3)*0.2)
    % % %
    caxis([0.05 1])
    % % %     % configurations to have the graphic ploted in the "normal" direction.
    axH = gca; %ax = gca returns the handle to the current axes for the
    % % %     % current figure. If an axes does not exist, then gca creates an axes and
    % % %     % returns its handle. You can use the axes handle to query and modify axes
    % % %     % properties.
    set(axH,'YDir','normal');
    % % %
    % % %     %%%%%%%%%%%%%%%removing number in the axis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %
    % %                  set(gca,'YTick',[])
    % %                 set(gca,'XTick',[])
    % % %
    % % %
    % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %increase dRindx
     rDindxN=rDindxN-1;
end

% % %
% % %include color bar
ax=gca;
pos=get(gca,'pos');
set(gca,'pos',[pos(1) pos(2) pos(3) pos(4)*0.95]);
pos=get(gca,'pos');
hc=colorbar('location','northoutside','position',[pos(1) pos(2)+pos(4) pos(3) 0.03]);
set(hc,'xaxisloc','top');
%
% %save figure
%
% %Save figure
figH = gcf;
set(figH,'Name',lRStr{lRindx});
resDir=[saveDir,filesep, dRT{1},filesep,rDtarget{1},aPtarget{1},lRtarget{1 }];
mkdir(resDir)
outFile =[resDir ,filesep,'pMatrixPlot'];
%
% %                     saveas(figH,outFile,'png');
saveas(figH,outFile,'fig');
%
fprintf('\nFigures saved in %s.\n',outFile);


%inclues a title for

[ax,h]=subtitle([rDtarget{1},aPtarget{1},lRtarget{1},dRT{1}]);
