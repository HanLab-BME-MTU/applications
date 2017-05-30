%script to plot the informations maximum p-value counts considering all 
%bootstrapping repetitions


%%
currDir ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/bootstrapping/results';
% the name until target
% title of the figure will consider the infos that are filled here for the
% name of the target
 rDtarget = {'rD10'};%,,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
 aPtarget = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
 lRtarget ={'lR0p3'};%,'lR0p3','lR0p4','lR0p5'};

   % values of rD, aP and lR of probe
    rDvals = [4;6;8;10;12;14;16];%20;40;60;80;100;140;160];
    aPvals = [0.2;0.3;0.4;0.5;0.6;0.7;0.8];
    lRvals = [0.1;0.2;0.3;0.4;0.5;0.6];%'lR0p01lR0p02';'lR0p02lR0p04';'lR0p03lR0p06
   lRStr = {'lR0p1';'lR0p2';'lR0p3';'lR0p4';'lR0p5';'lR0p6'};

for lRindx=1:length(lRvals)
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

%% plot

% Counts for the possible values of rD and aP in the diff. bootstrapping
% repetitions
temp= load([currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },filesep,'matrixMaxPvalueCount.mat']);
maxPvalueCount=temp.matrixMaxPValueCount;
maxSize=max(max(maxPvalueCount));
minSize=min(min(maxPvalueCount));
%the result will be 5 matrices, to have the global maximum value we need to
%put this values together as a line vector and then calculates the
%maximum/min
%limite of s.

limMaxS=max(reshape(maxSize,1,length(maxSize)));
limMinS=min(reshape(minSize,1,length(maxSize)));

 %load([currDir,filesep,'pMatrix',lRStr{lRindx},'.mat']);
 imagesc(aPvals,rDvals,maxPvalueCount(:,:,lRindx));
 
    h = colorbar;
 caxis([limMinS limMaxS])
    ylabel(h, 'Couts','FontSize',15) 
%   title(axH,['target:',rDtarget{1},filesep,aPtarget{1},filesep,lRtarget{1} ' and Probe:',lRStr{lRindx}],'FontSize',12);
 
 % configurations to have the graphic ploted in the "normal" direction.
        axH = gca; %ax = gca returns the handle to the current axes for the
% current figure. If an axes does not exist, then gca creates an axes and 
% returns its handle. You can use the axes handle to query and modify axes
% properties. 
         set(axH,'YDir','normal');
         set(axH,'FontSize',22);
% % label info        
        xlabel(axH,'Association probability','FontSize',22);        
        ylabel(axH,'Receptor Density','FontSize',22);
figH = gcf; 

      
 %% Save figure
         set(figH,'Name',lRStr{lRindx});
        outFile = [currDir,filesep,rDtarget{rDTIndx},aPtarget{aPTIndx},lRtarget{lRTIndx },filesep,'countsMaxPvalue',lRStr{lRindx},'_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');  
        
  end 
     end
      end
end