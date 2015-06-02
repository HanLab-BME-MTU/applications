%Script to calculate aggregation state from compound tracks. Compound
%tracks with appended aggregation state (in defaul and alternative format)
%will be saved in same directory as original compound tracks.
%
%Khuloud Jaqaman, May 2015

sourceRoot = '/project/biophysics/jaqaman_lab/interKinetics/ryirdaw/2014/09/090514/probeISruns';

%Define strings for directory hierarchy as needed
rDDir = {'rD4','rD6','rD8','rD10','rD12','rD14','rD16'};
aPDir = {'aP0p2','aP0p3','aP0p4','aP0p5','aP0p6','aP0p7','aP0p8'};
outDirNum = 1:5;
lRDir = {'lR0p1','lR0p2','lR0p3','lR0p4','lR0p5'};

%define intensity mean and stadnard deviation. Must match simulation input
%or experimentally-derived values
intensityQuantum = [1 0.3];

%Define number of label ratio
numLabelRatio = length(lRDir);

fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different runs
        for outDirIndx = 1 : length(outDirNum)
            
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
                    aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),...
                    filesep,lRDir{lRDirIndx}];
                
                %load compTracks
                tmp = load([currDir,filesep,'compTracks.mat']);
                compTracksIn = tmp.compTracks;
                clear tmp
                
                %get aggregation state
                [compTracksAggregState,segmentStat] = aggregStateFromCompTracks_new(...
                    compTracksIn,intensityQuantum);

                fprintf('\n   Out = %d, lR = %s ',outDirIndx,lRDir{lRDirIndx});
                
                %save tracks with aggregation state
                save([currDir,'/compTracksAggregState'],'compTracksAggregState','segmentStat','-v7.3');
                
                clear compTracksAggregState compTracksIn
                
                fprintf('... done.');
                
            end %for each labelRatio
            
        end %for each outDir
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
