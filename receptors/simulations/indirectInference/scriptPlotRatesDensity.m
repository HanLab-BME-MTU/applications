%Script to combine rate and density results from different equivalent
%"movies". Means and stds etc. will be saved in same directory as cell
%array of individual movie results
%
%Khuloud Jaqaman, June 2015

%Define strings for directory hierarchy as needed
sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20160712_Analysis/probe';

%Define strings for directory hierarchy as needed
rDDir = {'rD100'}; %,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}; %'dR0p5','dR2p0','dR5p0','aP0p2','aP0p5','aP0p8','dC0p05','dC0p2'};
lRDir = {'lR0p01','lR0p02','lR0p03','lR0p04','lR0p05','lR0p06'};


%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
            
            %read individual results
            tmp = load(fullfile(currDir,'ratesAndDensityComb_dt0p01_T10.mat'));
            ratesDensityComb(lRDirIndx,aPDirIndx,rDDirIndx) = tmp;
            
        end %for each labelRatio
        
    end %for each aP
    
end %for each rD

%plot them

plotColor = {'r','k','b','c','g'};
colorIndx = 0;

legendTitles = [];

MIN_MOV = 5;

figure, hold on

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
            colorIndx = colorIndx + 1;
            
            %get relevant output
            onRate = ratesDensityComb(lRDirIndx,aPDirIndx,rDDirIndx).rateOnPerClust;
            onRate(onRate(:,3) < MIN_MOV,:) = NaN;
            offRate = ratesDensityComb(lRDirIndx,aPDirIndx,rDDirIndx).rateOffPerClust;
            offRate(offRate(:,3) < MIN_MOV,:) = NaN;
            densityC = ratesDensityComb(lRDirIndx,aPDirIndx,rDDirIndx).densityPerClust;
            densityC(densityC(:,3) < MIN_MOV,:) = NaN;
            
            %subplot 1: On rate
            subplot(1,3,1)
            hold on
            plot(onRate(:,1),'o','LineStyle','none','Color',plotColor{colorIndx});
            myErrorbar(onRate(:,1),onRate(:,2));
            
            %subplot 2: Off rate
            subplot(1,3,2)
            hold on
            plot(offRate(:,1),'o','LineStyle','none','Color',plotColor{colorIndx});
            myErrorbar(offRate(:,1),offRate(:,2));
            
            %subplot 3: Density
            subplot(1,3,3)
            hold on
            plot(densityC(:,1),'o','LineStyle','none','Color',plotColor{colorIndx});
            myErrorbar(densityC(:,1),densityC(:,2));
            
            legendTitles{end+1} = [rDDir{rDDirIndx} ' ' aPDir{aPDirIndx} ' ' lRDir{lRDirIndx}];
           
        end %for each labelRatio
        
    end %for each aP
    
end %for each rD

legend(legendTitles)

clear

