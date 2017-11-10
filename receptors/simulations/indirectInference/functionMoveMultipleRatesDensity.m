function functionMoveMultipleRatesDensity(sourceRoot,destinationRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)

%Function to copy multiple rate and density results files to new directory for further
%compilation and analysis. New directory is specified by "destinationRoot."
%Hierarchy in destination directory follows that in source directory but
%without the "out" layer, which will be now incorporated in the results
%file name.
%
%INPUT
%  sourceRoot : path for the receptorInfoLabeled.
%  
%  destinationRoot: path where the rates and densities will be copied. 
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
%OUTPUT
% There is no output.
fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different labeling ratios
        for lRDirIndx = 1 : length(lRDir)
               for dRDirIndx = 1 : length(dRDir)
            
            %create destination directory
            destDir = [destinationRoot,filesep,rDDir{rDDirIndx},filesep,...
                dRDir{dRDirIndx},filesep,aPDir{aPDirIndx},filesep,lRDir{lRDirIndx},filesep,'ind'];
                 mkdir(destDir);
            
            %iterate through the different runs
            
          
            
            
            for outDirIndx = 1 : length(outDirNum)
                
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
               dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),filesep,lRDir{lRDirIndx}];
                
                %copy rates and densities file to new directory
                %append file name with number indicating movie #
                copyfile(fullfile(currDir,'ratesAndDensity_dt0p1_T10.mat'),...
                    fullfile(destDir,['ratesAndDensity_dt0p1_T10_' int2str(outDirNum(outDirIndx)) '.mat']));
                
            end %for each labelRatio
             
             end % for each dissociation rate
            
        end %for each outDir
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
