function functionCalcMultipleAggregState(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir,intensityQuantum)

%function to calculate multiple aggregation state from compound tracks. Compound
%tracks with appended aggregation state (in default and alternative formats)
%will be saved in same directory as original compound tracks.
%
%SYNOPSIS functionExtractMultipleCompTracksdR(sourceRoot,rDDir,aPDir,dRDir,outDirNum,lRDir)
%
%INPUT
%  sourceRoot : path for the receptorInfoLabeled.
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
%  intensityInfo: Row vector with unit intensity mean and standard
%                      deviation (e.g. the intensity of a single
%                      fluorophore labeling a single receptor). 
%
% OUTPUT
% There is no output, the results will be saved in the sourceRoot directory
% as compTracksAggregState.m
%
% Luciana de Oliveira, April 2017



fprintf('\n===============================================================');

%The top level directory is that of receptor density
for rDDirIndx = 1 : length(rDDir)
    
    tic    
    %Iterate through association probability values per density
    for aPDirIndx = 1 : length(aPDir)
        
        for dRDirIndx = 1 : length(dRDir)
        fprintf('\nProcessing rD = %s, aP = %s ',rDDir{rDDirIndx},aPDir{aPDirIndx});
        
        %iterate through the different runs
        for outDirIndx = 1 : length(outDirNum)
            
            %iterate through the different labeling ratios
            for lRDirIndx = 1 : length(lRDir)
                
                %name of current directory
                currDir = [sourceRoot,filesep,rDDir{rDDirIndx},filesep,...
               dRDir{dRDirIndx},filesep, aPDir{aPDirIndx},filesep,'out',int2str(outDirNum(outDirIndx)),filesep,lRDir{lRDirIndx}];
            
                
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
        
        end %for each dRDir
        
    end %for each aP
    
    elapsedTime = toc;
    fprintf('\nElapsed time for aP = %s is %g seconds.\n',aPDir{aPDirIndx},elapsedTime);
    
end %for each rD

fprintf('\n\nAll done.\n');

clear
