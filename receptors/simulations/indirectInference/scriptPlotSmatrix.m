%script to plot s matrix: it loads the file with 3D smatrix, calculates the
%maximum and minimum value of S to estabilish the limit values of color
%bar. It generates 5 different figures corresponding to the label ratio and
% saves them in the current file.


% directory with s matrix
currDir = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160817/results/S_PMatrix_sT25_dT0p1/target';% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
     rDtarget = {'rD4'};%;'rD6';'rD8';'rD10';'rD12';'rD14';'rD16'}; 
    aPtarget = {'aP0p5'};
    lRtarget = {'lR0p2'};
 
    
   % values of rD, aP and lR of probe
    rDvals = [4;6;8;10;12;14;16];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRStr = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5'};
 
  %making figures for different lR
  for lRindx=1:length(lRStr)
%load s matrix
temp= load([currDir,rDtarget{1},aPtarget{1},lRtarget{1},filesep,'sMatrix.mat']);
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

limMaxS=max(reshape(maxSize,1,5));
limMinS=min(reshape(minSize,1,5));

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
        outFile = [currDir,rDtarget{1},aPtarget{1},lRtarget{1},filesep,'sMatrix_',lRStr{lRindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
        
  end % make figures for different lR    
    