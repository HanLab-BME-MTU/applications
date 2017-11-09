%Script to combine rate and density results from different equivalent
%"movies". Means and stds etc. will be saved in same directory as cell
%array of individual movie results
%
%Khuloud Jaqaman, June 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161117/analysis';

%Define strings for directory hierarchy as needed
rDDir = {'rD4'};%,'rD8','rD10','rD12','rD14','rD16'};
aPDir = {'aP0p5'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
outDirNum = 1:10;
lRDir ={'lR0p2','lR0p4'};%,'lR0p3','lR0p4','lR0p5'};
BootNumber=1:100;
%','lR0p3','lR0p4','lR0p5','lR0p6','lR1p0'};
fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
           for  bootIndx=1:length(BootNumber)
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx},filesep,'ind',filesep,'bootstrapping',int2str(BootNumber(bootIndx))];
            %read individual results
            tmp = load(fullfile(currDir,'ratesAndDensityInd_dt0p1_T10.mat'));
            
            %call function to combine results
            [rateOnPerClust,rateOffPerClust,densityPerClust,paramVarCovMat,paramMatrix] = ...
                combineClusterRatesAndDensityNew(tmp.ratesDensityPerMovie);
            
            %make a Directory
            combDirBoot=[sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx},filesep,'bootstrapping',int2str(BootNumber(bootIndx))];
            mkdir(combDirBoot)
            %save combined results
            
          
            save([combDirBoot,'/ratesAndDensityComb_dt0p1_T10'],'rateOnPerClust',...
                'rateOffPerClust','densityPerClust','paramVarCovMat','paramMatrix','-v7.3');
            
           end %for each labelRatio
        end %for each bootstrapping number
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear

