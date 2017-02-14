%script to plot s matrix: it loads the file with 3D smatrix, calculates the
%maximum and minimum value of S to estabilish the limit values of color
%bar. It generates 5 different figures corresponding to the label ratio and
% saves them in the current file.


currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170126/target/results';

% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
 rDtarget = {'rD10'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p1'};%,'lR0p3','lR0p4','lR0p5'};

   % values of rD, aP and lR of probe
    rDvals = [4;6;8;10;12;14;16];%20;40;60;80;100;140;160];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRStr = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'};%'lR0p01lR0p02';'lR0p02lR0p04';'lR0p03lR0p06'
    
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

temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx},filesep,'sMatrix.mat']);
sMatrix=temp.sMatrix;
sMatrixLog=log10(sMatrix);

%defining limits for the color bar

%first calculate the maximum/min value for each receptor label

maxSize=max(max(sMatrixLog));
minSize=min(min(sMatrixLog));
%the result will be 5 matrices, to have the global maximum value we need to
%put this values together as a line vector and then calculates the
%maximum/min
%limite of s.

limMaxS=max(reshape(maxSize,1,length(maxSize)));
limMinS=min(reshape(minSize,1,length(maxSize)));

%plot the figure

 imagesc(aPvals,rDvals,sMatrixLog(:,:,lRindx));
 colorbar
 % configurations to have the graphic ploted in the "normal" direction.
        axH = gca; %ax = gca returns the handle to the current axes for the
% current figure. If an axes does not exist, then gca creates an axes and 
% returns its handle. You can use the axes handle to query and modify axes
% properties. 
        set(axH,'YDir','normal');
        set(axH,'FontSize',12);
% label info        
        xlabel(axH,'Association probability','FontSize',12);        
        ylabel(axH,'Receptor Density','FontSize',12);
       h = colorbar;
      caxis([limMinS limMaxS]) %limits for the colorbar
       ylabel(h, 'log10(Mahalanobis distance)') 
%title        
        title(axH,['target:',rDtarget{1},filesep,aPtarget{1},filesep,lRtarget{1} ' and Probe:',lRStr{lRindx}],'FontSize',12);

%Save figure
        figH = gcf;
        set(figH,'Name',lRStr{lRindx});
        outFile = [currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx},filesep,'sMatrix_',lRStr{lRindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
        
  end % make figures for different lR    
     end
     end
 end