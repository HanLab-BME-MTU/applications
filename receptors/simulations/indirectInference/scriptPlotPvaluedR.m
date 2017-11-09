%Script to plot p values: it takes the pmatrix and plot in a 2D graphic
%where y is the receptor density and x is the association probability.
% It generates 5 different figures corresponding to the label ratio and
% saves them in the current file.

% directory with pValue matrix
currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170508/results/rD4dR0p5aP0p5lR0p2';
 % name of the target
rDtarget = {'rD4'};%,'rD8','rD10','rD12','rD14','rD16'};
aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
dRT={'dR0p75'};%{'dR1'};%{'dR0p25','dR0p5','dR0p75','dR1','dR1p25','dR1p5','dR1p75','dR2'};%,'lR0p3','lR0p4','lR0p5'};
lRtarget={'lR0p2'};
% values of rD, aP and lR of probe
dRvals = 0.25:0.25:2;
rDvals=2:2:16;
aPvals = 0.2:0.1:0.8;
lRStr = {'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5','lR0p6'};%,'lR0p2','lR0p3','lR0p4','lR0p5','lR0p6'};
dRdir={'dR0p75'};
%{'lR0p01lR0p02lR0p03lR0p04'                                         % lRProbe = {'lR0p01','lR0p02','lR0p03','lR0p04';
%              'lR0p02lR0p04lR0p06lR0p08'
% 'lR0p03lR0p06lR0p09lR0p12'
% 'lR0p04lR0p08lR0p12lR0p16'
% 'lR0p06lR0p12lR0p18lR0p24'     }  ;                            %            'lR0p02','lR0p04','lR0p06','lR0p08';
%               'lR0p03lR0p06lR0p09lR0p12'                               %            'lR0p03','lR0p06','lR0p09','lR0p12';
%                'lR0p04lR0p08lR0p12lR0p16'                            %            'lR0p04','lR0p08','lR0p12','lR0p16';
%                  'lR0p06lR0p12lR0p18lR0p24'};

%figures

for lRindx=1:length(lRStr)
%     for lRTIndx = 1 : length(lRtarget)
        %
%         for dRDirIndx = 1 : length(dRT)
            %     tic
            %     %Iterate through association probability values per density
%             for aPTIndx = 1 : length(aPtarget)
                %
                %
                %         %iterate through the different labeling ratios
%                 for rDTIndx = 1 : length(rDvals)
                    
                    %load pvalue matrix
                    %  temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },filesep,'pMatrix_lR0p06lR0p12lR0p18lR0p24.mat']);
%                     temp= load([currDir,filesep, dRDir{dRDirIndx},filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },filesep,'pMatrix.mat']);
%                      temp= load([currDir,filesep,dRT{1},'Target',filesep,dRdir{1},filesep,rDtarget{1},aPtarget{1},lRtarget{1 },filesep,'pMatrix.mat']);
                    load([currDir,filesep,'pMatrix.mat']);
%                     pMatrix=temp.pMatrix;
                    pMatrix(pMatrix<0.05)=0.05;
                    %here it is using the colormap that I defined
                    %plot the figure
                    ax = gca;
                    load('MyColormaps','mycmap')
                    colormap(ax,mycmap)
                    imagesc(aPvals,dRvals,pMatrix(:,:,6,lRindx));
                    %to create the grid
                    line(repmat([0 2],9,1).',repmat(0.125:0.25:2.125,2,1),'Color',ones(1,3)*0.2)
                    line(repmat(0.15:0.1:0.85,2,1),repmat([0; 2.5],1,8),'Color',ones(1,3)*0.2)
                    
                    
                    %plot the figure
                    %  imagesc(aPvals,rDvals,pMatrix);
                    
                    
                    
                    colorbar
                    
                    % configurations to have the graphic ploted in the "normal" direction.
                    axH = gca; %ax = gca returns the handle to the current axes for the
                    % current figure. If an axes does not exist, then gca creates an axes and
                    % returns its handle. You can use the axes handle to query and modify axes
                    % properties.
                    set(axH,'YDir','normal');
                    set(axH,'FontSize',25);
                    % label info
                    xlabel(axH,'Association probability','FontSize',25);
                    ylabel(axH,'Dissociation Rate','FontSize',25);
                    caxis([0.05 1]) %limits for the colorbar
                    h = colorbar;
                    ylabel(h, 'p-value','FontSize',25)
                    
                    %title        rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx}
                    %title
                    %         title(axH,['target:',rDtarget{rDTIndx},filesep,aPtarget{aPTIndx},filesep,lRtarget{lRTIndx } ' and Probe:',lRStr{lRindx}],'FontSize',12);
                    
                    %Save figure
                     figH = gcf;
                    set(figH,'Name',lRStr{lRindx});
%                     outFile = [currDir,filesep,dRT{1},'Target',filesep,dRdir{1},filesep,rDtarget{1},aPtarget{1},lRtarget{1 },filesep,'pMatrix_',lRStr{lRindx},'_plot'];
outFile = [currDir,filesep,'pMatrix_',lRStr{lRindx},'_plot'];
                   saveas(figH,outFile,'png');
                    saveas(figH,outFile,'fig');
%                     
                    fprintf('\nFigures saved in %s.\n',outFile);
                end
%              end
%         end
%     end
% end