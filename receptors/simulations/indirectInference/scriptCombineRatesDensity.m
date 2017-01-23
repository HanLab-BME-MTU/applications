%Script to combine rate and density results from different equivalent
%"movies". Means and stds etc. will be saved in same directory as cell
%array of individual movie results
%
%Khuloud Jaqaman, June 2015

sourceRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170112/diffLabelRatio/analysis/probe';
%Define strings for directory hierarchy as needed
rDDir = {'rD120'};%,'rD8','rD10','rD12','rD14','rD16'};
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};%,'aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
outDirNum =1:30;
lRDir = {'lR0p14';'lR0p16';'lR0p18';'lR0p20';'lR0p22';'lR0p24'};%,'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5';

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

