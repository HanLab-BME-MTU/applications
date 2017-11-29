% This script calculates and plot the association rate from the receptor density from the simulations in function of
% receptor density

% sourceRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170316/superResol/staticData/analysis/probeIS_sT25_dT0p1';
sourceRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/varyDissRate/analysis/probeIS_sT25_dT0p1';
saveRoot='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170407/superRes/plotOnRate/varyReceptorDensity/check';

%Define strings for directory hierarchy as needed
rDDir ={'rD16'};%,'rD20','rD40','rD60','rD80','rD100','rD120','rD140''rD8',
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%,'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}
outDirNum =1:10;
lRDir = {'lR1'};
dR=1.75;
dRDir={'dR1p75'};
%{'lR0p01';'lR0p02';'lR0p03';'lR0p04';'lR0p05';'lR0p06';'l%,'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5';

% save space for the onRate matrix
% onRate=zeros(1,length(aPDir));
onRate=zeros(1,length(aPDir));
legendTitles = [];
%  for dRindex=1: length(dRDir) 
for rDDirIndx = 1 : length(rDDir)
    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        for lRIndex=1:size(lRDir)
          
%            
                 currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
               dRDir{1},filesep,  aPDir{aPDirIndx},filesep,lRDir{lRIndex}];
%   currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
%               aPDir{aPDirIndx},filesep,lRDir{lRIndex}];
            
%             
            tmp =  load([fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat')]);
            paramMatrix=tmp.paramMatrix;
            
            %calculate the value of association rate
            
%             onRate(aPDirIndx)=dR(dRindex)*mean(paramMatrix(2,:))/(mean(paramMatrix(1,:)))^2;
             onRate(aPDirIndx)=dR*mean(paramMatrix(2,:))/(mean(paramMatrix(1,:)))^2;
%            end
        end
    end
    
    %plot
     plot(0.2:0.1:0.8,onRate,'o')
%       plot(0.2:0.1:0.8,onRate)
     xlim([0.2 0.8])
    set(gca,'fontsize',25)
    
     xlabel('Association probability','FontSize',25);        
         ylabel('On rate','FontSize',25);
%        legendTitles{end+1}=[rDDir{rDDirIndx}];
 

legendTitles{end+1}=rDDir{rDDirIndx};
    hold all
  
    
    
    h_legend= legend(legendTitles,'Location','northwest');
   set(h_legend,'FontSize',10);
end

%Save figure
%         figH = gcf;
%         outFile = [saveRoot,filesep,'onRateVaryDr_',dRDir{1}];
%         saveas(figH,outFile,'png');
%         saveas(figH,outFile,'fig');
%         
%         fprintf('\nFigures saved in %s.\n',outFile);