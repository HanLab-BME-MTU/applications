%Script to combine rate and density results from different equivalent
%"movies". Combined results are stored as a cell array and saved in same
%directory next to individual results.
%
%Khuloud Jaqaman, June 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/kjaqaman/150608_AnalysisLR1_dt0p01_T10/probeISruns';

%Define strings for directory hierarchy as needed
rDDir = {'rD10'}; %,'rD40','rD60','rD80','rD100','rD120','rD140','rD160'};
aPDir = {'dR0p5','dR2p0','dR5p0','aP0p2','aP0p5','aP0p8','dC0p05','dC0p2'};
outDirNum = 1:10;
lRDir = {'lR1p0'}; %,'lR0p2','lR0p3','lR0p4','lR0p5','lR0p6','lR1p0'};

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
            
            %define output cell array
            ratesDensityPerMovie = cell(length(outDirNum),1);
            
            %iterate through the different runs
            for outDirIndx = 1 : length(outDirNum)
                
                %read results
                tmp = load(fullfile(currDir,'ind',['ratesAndDensity_dt0p01_T10_' int2str(outDirNum(outDirIndx)) '.mat']));
                
                %store in cell array
                ratesDensityPerMovie{outDirIndx} = tmp;
                
            end %for each outDir
            
            %save cell array
            save([currDir,'/ratesAndDensityInd_dt0p01_T10'],'ratesDensityPerMovie','-v7.3');
            
        end %for each labelRatio
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
                
