function functionCombineMultopleRatesDensityStaticDynamicdR(destinationRoot,rDDir,aPDir,dRDir,lRDir)
 
%function to combine rate and density results from different equivalent
%"movies". Means and stds etc. will be saved in same directory as cell
%array of individual movie results
%INPUT
%  destinationRoot: path to ratesAndDensity_dt0p1_T10. 
%
%  rDDir: multiple receptor density directories that will be analyzed. It
%  is a cell in the form: {'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'}
%
%  aPDir      : association probability directories. A cell in the form: {'aP0p2','aP0p3','aP0p4'}
%
%  dRDir      : dissociation rate directories. A cell in the form: {'aP0p2','aP0p3','aP0p4'}
%
%  outDirNum  : list with the simulation number. 1:10 for the low density and
%  1:30 for the high density.
%
%  lRDir      : labed receptor directory. A cell in the form: {'lR0p2';'lR0p3';'lR0p4'}
%
% OUTPUT
% There is no output, the results will be saved in the sourceRoot directory
% as ratesAndDensityComb_dt0p1_T10.m
fprintf('\n===============================================================');



%%%%%%%%%%%%%%%%%%%%%%%%%
% here is the definition if it is static (0) or dynamic (1) data
systemState=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    for dRDirIndx = 1 : length(dRDir)
        tic
        %Iterate through association probability values per density
        for aPDirIndx = 1 : length(aPDir)
            
            fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
            
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                
                %name of current directory
                currDir = [destinationRoot,filesep,rDDir{rDDirIndx},filesep,...
                    dRDir{dRDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,lRDir{lRDirIndx}];
                
                %read individual results
                
                
                
                tmp = load(fullfile(currDir,'ratesAndDensityInd_dt0p1_T10.mat'));
                for iboot=1: bootRepSize
                ratesDensityPerMovie=tmp.ratesDensityPerMovie;
                %call function to combine results
                [rateOnPerClust{bootRepSize},rateOffPerClust{bootRepSize},densityPerClust{bootRepSize},paramVarCovMat{bootRepSize},paramMatrix{bootRepSize}] = ...
                    combineClusterRatesAndDensityStaticDynamic(ratesDensityPerMovie{bootRepSize},systemState);
                %save combined results
               
                end
                save([currDir,'/ratesAndDensityComb_dt0p1_T10'],'rateOnPerClust',...
                    'rateOffPerClust','densityPerClust','paramVarCovMat','paramMatrix','-v7.3');
                
                
            end %for each labelRatio
            
        end %for each aP
        
        elapsedTime = toc;
        fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
        
    end %for each rD
end

fprintf('\n\nAll done.\n');

clear

