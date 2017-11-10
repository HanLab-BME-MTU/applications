%Script to plot the receptor fraction in function of frames

%Define strings for directory hierarchy as needed
sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2017/01/20170112/diffLabelRatio/target';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170331/equilibriumState/out1';
%Define strings for directory hierarchy as needed
rDDir = {'rD100'}; %,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
aPDir = {'aP0p5'}; %'dR0p5','dR2p0','dR5p0','aP0p2','aP0p5','aP0p8','dC0p05','dC0p2'};
outDirNum=1;
lRDir = {'lR0p14','lR0p16','lR0p18','lR0p20','lR0p22','lR0p24'};
dRDir={'dR1'};%,'dR1p5','dR1p75' 'dR0p5','dR2'

% load the clustStatus-clustFrac

%name of current directory
for dRindex=1:length(dRDir);
    for lRIndex=1:length(lRDir);
    
currDir = [sourceRoot,filesep,rDDir{1},filesep,...
    aPDir{1},filesep,'out',int2str(outDirNum(1)),filesep,lRDir{lRIndex}];

%read individual results
tmp = load(fullfile(currDir,'ratesAndDensity_dt0p1_T10.mat'));
clustStats=tmp.clustStats;
clustFrac=clustStats.clusterFrac;

% calculate the maximum cluster size
maxClustSize=size(clustFrac,1);

%plot

figure, hold on
plotColor = {'r','c','k','b','g'};
colorIndx = 0;
MIN_MOV = 5;
legendTitles = [];
for clustSizeIndex=1:maxClustSize
    if length( clustFrac(clustSizeIndex,:))>MIN_MOV
     colorIndx = colorIndx + 1;
 plot(clustFrac(clustSizeIndex,:),'Color',plotColor{colorIndx}); 
 axis([0 101 0 1])
    end
legendTitles{end+1}=['cluster size ', int2str(clustSizeIndex)];
xlabel('frames','FontSize',15);        
  
ylabel('Receptor density fraction','FontSize',15);
end
legend(legendTitles)

   
%Save figure
        figH = gcf;
        outFile = [saveRoot,filesep,rDDir{1},aPDir{1},lRDir{lRIndex},dRDir{dRindex},'clustFrac_plot'];
        saveas(figH,outFile,'png');
        saveas(figH,outFile,'fig');
        
        fprintf('\nFigures saved in %s.\n',outFile);
    end
end
