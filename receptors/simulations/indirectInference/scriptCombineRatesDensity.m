%Script to combine rate and density results from different equivalent
%"movies". Means and stds etc. will be saved in same directory as cell
%array of individual movie results
%
%Khuloud Jaqaman, June 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20161004/analysis';

%Define strings for directory hierarchy as needed
rDDir = {'rD20'};%,'rD60','rD80','rD120','rD140','rD160'}; aPDir = {'aP0p5'};
lRDir = {'lR0p03'};
 %,'lR0p02','lR0p03','lR0p04','lR0p05','lR0p06'};

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
            %name of current directory
            currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
            
            %read individual results
            tmp = load(fullfile(currDir,'ratesAndDensityInd_dt0p1_T10.mat'));
            
            %call function to combine results
            [rateOnPerClust,rateOffPerClust,densityPerClust,paramVarCovMat] = ...
                combineClusterRatesAndDensity(tmp.ratesDensityPerMovie);
            
            %save combined results
            save([currDir,'/ratesAndDensityComb_dt0p1_T10'],'rateOnPerClust',...
                'rateOffPerClust','densityPerClust','paramVarCovMat','-v7.3');
            
        end %for each labelRatio
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear

