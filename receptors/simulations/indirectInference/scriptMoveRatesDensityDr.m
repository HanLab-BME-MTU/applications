%Function to copy multiple rate and density results files to new directory for further
%compilation and analysis. New directory is specified by "destinationRoot."
%Hierarchy in destination directory follows that in source directory but
%without the "out" layer, which will be now incorporated in the results
%file name.
%
%Khuloud Jaqaman, June 2015

sourceRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/20170424/targetISruns';
destinationRoot ='/project/biophysics/jaqaman_lab/interKinetics/ldeoliveira/2016/12/20161201/analysis/targetIS_sT25_dT0p1';
%Define strings for directory hierarchy as needed
 rDDir = {'rD4'};%,'rD2','rD4','rD6','rD8','rD10','rD12','rD14','rD16'
 aPDir = {'aP0p5'};%'aP0p6','aP0p7','aP0p8'
 dRDir={'dR1'};
outDirNum =1:10;
lRDir = {'lR0p4'};



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
