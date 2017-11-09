%script to plot s matrix: it loads the file with 3D smatrix, calculates the
%maximum and minimum value of S to estabilish the limit values of color
%bar. It generates 5 different figures corresponding to the label ratio and
% saves them in the current file.


currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161202/lowRDProbes-highRDTargets_comparison/results/pairComparison/target_sT25_dT0p1';
saveDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161202/lowRDProbes-highRDTargets_comparison/results/pairComparison/target_sT25_dT0p1';
% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
 rDtarget = {'rD10'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p2lR0p4'};%,'lR0p3','lR0p4','lR0p5'};

   % values of rD, aP and lR of probe
    rDvals = [20;40;60;80;100;120;140;160];%20;40;60;80;100;120;140;160;
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRStr = {'lR0p01lR0p02';'lR0p02lR0p04';'lR0p03lR0p06'};%{'lR0p01','lR0p02','lR0p03','lR0p04','lR0p05','lR0p06'};
   %  {'lR0p01lR0p02';'lR0p02lR0p04';'lR0p03lR0p06'} 
   
 %figures   

 for lRindx=1:length(lRStr)
     
     for lRTIndx = 1 : length(lRtarget)
%     
%     
%     tic
     %Iterate through association probability values per density
     for aPTIndx = 1 : length(aPtarget)
%         
%                
%         %iterate through the different labeling ratios
         for rDTIndx = 1 : length(rDtarget)

    %load pvalue matrix

temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx},filesep,'sMatrix',lRStr{lRindx},'.mat']);
sMatrix=temp.sMatrix;
 sMatrixLog=log10(sMatrix);

%defining limits for the color bar

%first calculate the maximum/min value for each receptor label

maxSize(lRindx)=max(max(sMatrixLog));
minSize(lRindx)=min(min(sMatrixLog));
limMaxS=max(maxSize);
limMinS=min(minSize);
         end
     end
     end
 end
%the result will be 5 matrices, to have the global maximum value we need to
%put this values together as a line vector and then calculates the
%maximum/min
%limite of s.

% limMaxS=max(reshape(maxSize,1,length(maxSize)));
% limMinS=min(reshape(minSize,1,length(maxSize)));

%plot the figure
for lRindx=1:length(lRStr)
     
     for lRTIndx = 1 : length(lRtarget)
%     
%     
%     tic
     %Iterate through association probability values per density
     for aPTIndx = 1 : length(aPtarget)
%         
%                
%         %iterate through the different labeling ratios
         for rDTIndx = 1 : length(rDtarget)
             
             temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx},filesep,'sMatrix',lRStr{lRindx},'.mat']);
sMatrix=temp.sMatrix;
 sMatrixLog=log10(sMatrix);
 imagesc(aPvals,rDvals,sMatrixLog);
 colorbar
 % configurations to have the graphic ploted in the "normal" direction.
        axH = gca; %ax = gca returns the handle to the current axes for the
% current figure. If an axes does not exist, then gca creates an axes and 
% returns its handle. You can use the axes handle to query and modify axes
% properties. 
        set(axH,'YDir','normal');
        set(axH,'FontSize',15);
% label info        
        xlabel(axH,'Association probability','FontSize',16);        
        ylabel(axH,'Receptor Density','FontSize',16);
       h = colorbar;
       caxis([limMinS limMaxS]) %limits for the colorbar
       ylabel(h, 'log10(Mahalanobis distance)','FontSize',15) 
%title        
        title(axH,['target:',rDtarget{1},filesep,aPtarget{1},filesep,lRtarget{1} ' and Probe:',lRStr{lRindx}],'FontSize',15);

%Save figure
        figH = gcf;
        set(figH,'Name',lRStr{lRindx});
        outFile = [saveDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx},filesep,'sMatrix_',lRStr{lRindx},'_plot_1'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
        
  end % make figures for different lR    
     end
     end
 end