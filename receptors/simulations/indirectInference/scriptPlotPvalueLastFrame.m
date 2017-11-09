%Script to plot p values: it takes the pmatrix and plot in a 2D graphic
%where y is the receptor density and x is the association probability. 
% It generates 5 different figures corresponding to the label ratio and
% saves them in the current file.

% directory with pValue matrix
currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170207/results';

% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
 rDtarget = {'rD100'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p03'};%,'lR0p3','lR0p4','lR0p5'};

   % values of rD, aP and lR of probe
    rDvals = [20;40;60;80;100;120;140];%20;40;60;80;100;140;160];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRStr = {'lR0p03'};%{'lR0p08';'lR0p09';'lR0p1';'lR0p12';'lR0p13';'lR0p06'}
    
 %figures   
 
  for lRindx=1:length(lRStr)
      for lRTIndx = 1 : length(lRtarget)
%     
%     
%     tic
%     %Iterate through association probability values per density
     for aPTIndx = 1 : length(aPtarget)
%         
%                
%         %iterate through the different labeling ratios
  for rDTIndx = 1 : length(rDtarget)

    %load pvalue matrix
 temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },filesep,'pMatrix',lRStr{lRindx},'.mat']);
 %load([currDir,filesep,'pMatrix',lRStr{lRindx},'.mat']);
 pMatrix=temp.pMatrix;
% pMatrix(pMatrix<0.05)=0.05;
%plot the figure

 imagesc(aPvals,rDvals,pMatrix);
 colorbar
 
 % configurations to have the graphic ploted in the "normal" direction.
        axH = gca; %ax = gca returns the handle to the current axes for the
% current figure. If an axes does not exist, then gca creates an axes and 
% returns its handle. You can use the axes handle to query and modify axes
% properties. 
        set(axH,'YDir','normal');
        set(axH,'FontSize',22);
% label info        
        xlabel(axH,'Association probability','FontSize',22);        
         ylabel(axH,'Receptor Density','FontSize',22);
        caxis([0 1]) %limits for the colorbar
       h = colorbar;
       ylabel(h, 'p-value','FontSize',22) 
       
%title        rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx}
%title        
%         title(axH,['target:',rDtarget{rDTIndx},filesep,aPtarget{aPTIndx},filesep,lRtarget{lRTIndx } ' and Probe:',lRStr{lRindx}],'FontSize',12);
        
%Save figure
        figH = gcf;
        set(figH,'Name',lRStr{lRindx});
        outFile = [currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },filesep,'pMatrix_',lRStr{lRindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
  end      
     end 
      end
  end