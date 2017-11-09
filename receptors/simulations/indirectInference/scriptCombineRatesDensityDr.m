%Function to combine rate and density results from different equivalent
%"movies". Means and stds etc. will be saved in same directory as cell
%array of individual movie results
%

sourceRoot= '/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161201/analysis/targetIS_sT25_dT0p1';
%Define strings for directory hierarchy as needed
 rDDir = {'rD4'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
 aPDir = {'aP0p5'};%'aP0p6','aP0p7','aP0p8'
 dRDir={'dR1'};
outDirNum =1:10;
lRDir = {'lR0p4'};
%,'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5';
%,'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5';
fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
         
        for dRDirIndx = 1 : length(dRDir)
            
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
            
            %name of current directory
             currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
              dRDir{dRDirIndx},filesep,  aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
            
            
            %read individual results
            tmp = load(fullfile(currDir,'ratesAndDensityInd_dt0p1_T10.mat'));
            
           %call function to combine results
            [rateOnPerClust,rateOffPerClust,densityPerClust,paramVarCovMat,paramMatrix] = ...
                combineClusterRatesAndDensity(tmp.ratesDensityPerMovie);
            
            %save combined results
            save([currDir,'/ratesAndDensityComb_dt0p1_T10'],'rateOnPerClust',...
                'rateOffPerClust','densityPerClust','paramVarCovMat','paramMatrix','-v7.3');
        end
            
        end %for each labelRatio
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear

