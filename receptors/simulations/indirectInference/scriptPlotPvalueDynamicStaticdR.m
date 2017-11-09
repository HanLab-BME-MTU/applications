%Script to plot p values: it takes the pmatrix and plot in a 2D graphic
%where y is the receptor density and x is the association probability. 
% It generates 5 different figures corresponding to the label ratio and
% saves them in the current file.

% directory with pValue matrix
currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170508/results/combinationDynamicStatic/targetdR1p5';
% 
% 
% % name of the target
 rDtarget = {'rD12'};%'rD20','rD40','rD60','rD80','rD100','rD120','rD140','rD160'
 aPtarget = {'aP0p6'};%'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'
lRtarget = {'lR1'};%,'lR0p3','lR0p4','lR0p5'};
dRtarget={'dR1p5'};
dRProbe={'dR1p5'};


   % values of rD, aP and lR of probe
    rDvals = 2:2:16;
    aPvals = 0.2:0.1:0.8;
     lRStr = {'lR1'};
%      lRStr = {'lR0p1lR1'
%  'lR0p2lR1'
%  'lR0p3lR1'
%  'lR0p4lR1'
%  'lR0p5lR1'
%  'lR0p6lR1'};%{'lR0p01','lR0p02','lR0p03','lR0p04','lR0p05','lR0p06'}
    
 %figures   
 
  for lRindx=1:length(lRStr)
      
      for lRTIndx = 1 : length(lRtarget)
%     
     for dRIndex=1:length(dRProbe) 
%     tic
%     %Iterate through association probability values per density
     for aPTIndx = 1 : length(aPtarget)
%         
%                
%         %iterate through the different labeling ratios
  for rDTIndx = 1 : length(rDtarget)

    %load pvalue matrix
%  temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },dRProbe{dRIndex},filesep,'pMatrix_',lRStr{1},'.mat']);
% temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },dRProbe{dRIndex},filesep,'pMatrix_',lRStr{1},'.mat']);

temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },dRProbe{dRIndex},filesep,'pMatrix_',lRStr{1},'.mat']);

 %load([currDir,filesep,'pMatrix',lRStr{lRindx},'.mat']);
 %load([currDir,filesep,'pMatrix',lRStr{lRindx},'.mat']);
 pMatrix=temp.pMatrix;
pMatrix(pMatrix<0.05)=0.05;
 %plot the figure
                    ax = gca;
                    load('MyColormaps','mycmap')
                    colormap(ax,mycmap)
                    imagesc(aPvals,rDvals,pMatrix(:,:,lRindx));
                    %to create the grid
                    line(repmat([0 1],9,1).',repmat(1:2:17,2,1),'Color',ones(1,3)*0.2)
                    line(repmat(0.15:0.1:0.85,2,1),repmat([1; 17],1,8),'Color',ones(1,3)*0.2)
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
         ylabel(axH,'Receptor Density','FontSize',25);
        caxis([0.05 1]) %limits for the colorbar
       h = colorbar;
       ylabel(h, 'p-value','FontSize',25) 
       
%title        rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx}
%title        
%          title(axH,['target:',rDtarget{rDTIndx},filesep,aPtarget{aPTIndx},filesep,lRtarget{lRTIndx },filesep, dRtarget{1} ' and Probe:',lRStr{lRindx}, filesep, dRProbe{1}],'FontSize',12);
        
%Save figure
        figH = gcf;
        set(figH,'Name',lRStr{lRindx});
%         outFile = [currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },dRProbe{dRIndex},filesep,'pMatrix_',lRStr{lRindx},'_plot'];
outFile = [currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },dRProbe{dRIndex},filesep,'pMatrix_',lRStr{lRindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
  end
  fprintf('\nFigures saved in %s.\n',outFile);
  end      
     end 
      end
  end