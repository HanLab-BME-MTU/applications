%Script to plot target and probe from different folders using the outputs
%from scriptCombineRatesDensity

%Luciana de Oliveira July 2016

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/kjaqaman/150605_AnalysisAll_dt0p1_T10/probeISruns';
sourceRoot1 = '/project/biophysics/jaqaman_lab/interKinetics/kjaqaman/150605_AnalysisAll_dt0p1_T10/targetISruns';

%Define strings for directory hierarchy as needed
rDDir = {'rD6'}; %,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
aPDir = {'aP0p3'}; %,'aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'}; %'dR0p5','dR2p0','dR5p0','aP0p2','aP0p5','aP0p8','dC0p05','dC0p2'};
lRDir = {'lR0p2'}; %,'lR0p2','lR0p3','lR0p4','lR0p5','lR0p6','lR1p0'};

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
             currDir1 = [sourceRoot1,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
            %read individual results
            tmp = load(fullfile(currDir,'ratesAndDensityComb_dt0p1_T10.mat'));
            ratesDensityComb(lRDirIndx,aPDirIndx,rDDirIndx) = tmp;
            tmp1 = load(fullfile(currDir1,'ratesAndDensityComb_dt0p1_T10.mat'));
            ratesDensityComb1(lRDirIndx,aPDirIndx,rDDirIndx) = tmp1;
        end %for each labelRatio
        
    end %for each aP
    
end %for each rD

%plot them

plotColor = {'r','k','b','c','g','m'};
colorIndx = 0;
plotColor1 = {'k','b','k','c','g','m'};
colorIndx1 = 0;
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
            
            onRate1 = ratesDensityComb1(lRDirIndx,aPDirIndx,rDDirIndx).rateOnPerClust;
            onRate1(onRate1(:,3) < MIN_MOV,:) = NaN;
            offRate1 = ratesDensityComb1(lRDirIndx,aPDirIndx,rDDirIndx).rateOffPerClust;
            offRate1(offRate1(:,3) < MIN_MOV,:) = NaN;
            densityC1 = ratesDensityComb1(lRDirIndx,aPDirIndx,rDDirIndx).densityPerClust;
            densityC1(densityC1(:,3) < MIN_MOV,:) = NaN;
            %subplot 1: On rate
            subplot(1,3,1) %subplot(m,n,p) divides the current figure
%           into an m-by-n grid and creates an axes for a subplot in the 
%           position specified by p. MATLAB® numbers its subplots by row,
%           such that the first subplot is the first column of the first row,
%           the second subplot is the second column of the first row, and so on.
%           If the axes already exists, then the command subplot(m,n,p) makes
%           the subplot in position p the current axes.
            hold on
             plot(onRate(:,1),'o','LineStyle','none','Color',plotColor{colorIndx});
             myErrorbar(onRate(:,1),onRate(:,2));
            plot(onRate1(:,1),'o','LineStyle','none','Color',plotColor1{colorIndx});
            myErrorbar(onRate1(:,1),onRate1(:,2));
             axis([1.95 2.05 0 0.06])
            
           
            xlabel('Oligomer size')
            ylabel('Aggregation rate (per s per #/ \mu m^2)')
            %subplot 2: Off rate
            subplot(1,3,2)
            hold on
            plot(offRate(:,1),'o','LineStyle','none','Color',plotColor{colorIndx});
            myErrorbar(offRate(:,1),offRate(:,2));
            plot(offRate1(:,1),'o','LineStyle','none','Color',plotColor1{colorIndx});
            myErrorbar(offRate1(:,1),offRate1(:,2));
            axis([1.95 2.05 0 1])
            xlabel('Oligomer size')
            ylabel('Dissociation rate (per s per #/ \mu m^2)')
            
            
            %subplot 3: Density
            subplot(1,3,3)
            hold on
            plot(densityC(:,1),'o','LineStyle','none','Color',plotColor{colorIndx});
            myErrorbar(densityC(:,1),densityC(:,2));
            plot(densityC1(:,1),'o','LineStyle','none','Color',plotColor1{colorIndx});
            myErrorbar(densityC1(:,1),densityC1(:,2));
            xlabel('Oligomer size')
            ylabel('Density (#/ \mu m^2)')
            
            
            legendTitles{end+1} = [rDDir{rDDirIndx} ' ' aPDir{aPDirIndx} ' ' lRDir{lRDirIndx}];
            
            legend('probe:rD6 ap0p3 lR0p2 ' ,'target:rD6 ap0p3 lR0p2')
        end %for each labelRatio
        
    end %for each aP
    
end %for each rD



clear

