function collectParSimOutput(numSims,rootDir,outFileID)
%COLLECTPARSIMOUTPUT collects and saves simulation output values from a set 
%of parallel runs.
%
%   Values from the individual run will be accumulated in cells or arrays
%   of the same name with the string 'All' appended to the name.  Arrays
%   are used for the three rates calculated (association, association mono
%   and dissociation).  Also, recept2clustAssign and mismatchInfo are saved
%   by themselves while all other variables are saved in a single file name
%   valuesFor_<outFileID>.mat.  The first two combined with the rest of the
%   output variables can add up to over 2 GB. To be safe, use '-v7.3'
%   option when saving.
%
%   INPUT:
%       numSim:     interger number of simulations
%       rootDir:    string for root location of the output files.  For a given
%                   simulation i, the output files are expected to be in a 
%                   file named outi.mat in a folder named outi, under the
%                   root directory rootDir.
%       outFileID:  identifier string for output file with collected
%                    values to be written.
%
%   OUTPUT:
%       The following collected quantities are saved in mat files:
%       compTracksALTAll:       set of compTracksALT
%       mismatchInfoAll:        set of mismatchInfo
%       recept2clustAssignAll:  set of recept2clustAssign
%       valueFor_outFileID:     all other quantites saved together
%
%   Robel Yirdaw, 03/11/14
%       Modified, 05/20/14
%       Modified, 07/09/14
%       Modified, 08/26/14
%       Modified, 09/10/14
%

    %Flag to check for input arguments
    inputErr = 0;
    %All three required
    if (nargin ~= 3)
        fprintf('\nNumber of runs, root directory and output file id required.');
        inputErr = 1;
    end
    
    %Continue on if correct input arguments received
    if(~inputErr)
        
        %Flag for any error that might occur during file read/write and
        %processing of the simulation output files.
        valCollectionErr = 0;
        simIndx = NaN;
        
        try            
            %Change to root location of output folders
            cd(rootDir);
            
            %Allocate the necessary memory
            aggregStateAll = cell(numSims,1);
            aggregStateFromR2CAll = cell(numSims,1);
            clustHistoryMergedAll = cell(numSims,1);
            eventTableAll = cell(numSims,1);
            eventTable_monoAll = cell(numSims,1);
            mismatchInfoAll = cell(numSims,1);
            modelParamAll = cell(numSims,1);
            rateAssocPerClustAll = NaN(1000,numSims);
            rateDissocPerClustAll = NaN(1000,numSims);
            rateAssocMonoAll = NaN(1000,numSims);
            recept2clustAssignAll = cell(numSims,1);
            segmentStatAll = cell(numSims,1);
            simParamAll = cell(numSims,1);
            %070914
            elapsedTimeAll = NaN(numSims,1);
            assocStatsAll = cell(numSims,1);
            %082614
            collProbStatsAll = cell(numSims,1);
            %091014 - currently not collecting compTracksALT for indirect
            %inference work
            %compTracksALTAll = cell(numSims,1);
            
            %Will track the largest cluster size
            largestClust = 0;
            %Will track the largest cluster size with at least 10
            %associtiation and dissociation events - to be used when
            %calculating means of rates
            cutOffClust = inf;
            
            %Process all simulation outputs
            for simIndx=1:numSims
                fprintf('\nProcessing simulation #%d...',simIndx);
                
                currSimVals = load(['out',int2str(simIndx),'/out',int2str(simIndx),'.mat']);

                if (isfield(currSimVals,'aggregStateALL'))
                    aggregStateAll{simIndx} = currSimVals.aggregStateALL;
                end
                if (isfield(currSimVals,'aggregStateFromR2C'))
                    aggregStateFromR2CAll{simIndx} = currSimVals.aggregStateFromR2C;
                end
                clustHistoryMergedAll{simIndx} = currSimVals.clustHistoryMerged;
                eventTableAll{simIndx} = currSimVals.eventTable;
                eventTable_monoAll{simIndx} = currSimVals.eventTable_mono;                 
                if (isfield(currSimVals,'mismatchInfo'))
                    mismatchInfoAll{simIndx} = currSimVals.mismatchInfo;
                end
                modelParamAll{simIndx} = currSimVals.modelParam;

                currLargestClust = length(currSimVals.eventTable(:,1));
                rateAssocPerClustAll(1:currLargestClust,simIndx) = ...
                    currSimVals.rateAssocPerClust;
                rateDissocPerClustAll(1:currLargestClust,simIndx) = ...
                    currSimVals.rateDissocPerClust;
                rateAssocMonoAll(1:length(currSimVals.rateAssocMono),simIndx) = ...
                    currSimVals.rateAssocMono;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %052014 - the values in recept2clustAssign should not
                %exceed intmax('uint16') or 65535 so save as uint16.
                if (isfield(currSimVals,'recept2clustAssign'))
                    recept2clustAssignAll{simIndx} =...
                        uint16(currSimVals.recept2clustAssign);
                else
                    loadStruct = load(['out',int2str(simIndx),'/receptorInfoLabeled',int2str(simIndx),'.mat']);
                    recept2clustAssignAll{simIndx} =...
                        uint16(loadStruct.receptorInfoLabeled.recept2clustAssign);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if (isfield(currSimVals,'segmentStat'))
                    segmentStatAll{simIndx} = currSimVals.segmentStat;
                end
                simParamAll{simIndx} = currSimVals.simParam;
                
                %070914
                if (isfield(currSimVals,'elapsedTime'))
                    elapsedTimeAll(simIndx,1) = currSimVals.elapsedTime;
                end
                if (isfield(currSimVals,'assocStats'))
                    assocStatsAll{simIndx} = currSimVals.assocStats;
                end
                %082614
                if (isfield(currSimVals,'collProbStats'))
                    collProbStatsAll{simIndx} = currSimVals.collProbStats;
                end
                %091014
                %tempCompTracksALT = load(['out',int2str(simIndx),'/compTracksALT',int2str(simIndx),'.mat']);
                %compTracksALTAll{simIndx} = tempCompTracksALT.compTracksALT;

                %Update largest cluster size if necessary
                if (currLargestClust > largestClust)
                    largestClust = currLargestClust;
                end
                
                %Find the cut-off cluster and update the saved value if
                %necessary
                currCutOffClust = find(( ( currSimVals.eventTable(:,3) < 10) |...
                    (currSimVals.eventTable(:,6) < 10) ),1,'first');
                if (currCutOffClust < cutOffClust)
                    cutOffClust = currCutOffClust;
                end

                clear currLargestClust currCutOffClust currSimVals 
                
                fprintf('done.');
            end
            
        catch valCollectionExcep

            fprintf('\nException copying values (at simulation #%d).\n',simIndx);
            disp(valCollectionExcep.message);

            valCollectionErr = 1;
        end
        
        %If no error encountered during processing, write accumulated
        %variables to current folder (rootDir)
        if (~valCollectionErr)
            
            fprintf('\n\nWriting collected values to file...');
            
            %Trim rate matrices
            rateAssocPerClustAll(largestClust+1:end,:) = [];
            rateDissocPerClustAll(largestClust+1:end,:) = [];
            rateAssocMonoAll(largestClust:end,:) = [];
            
            %Write these two separately and clear afterwards
            save(['mismatchInfoAll_',outFileID,'.mat'],'mismatchInfoAll','-v7.3');
            save(['recept2clustAssignAll_',outFileID,'.mat'],'recept2clustAssignAll','-v7.3');
            %091014
            %save(['compTracksALTAll_',outFileID,'.mat'],'compTracksALTAll','-v7.3');
            
            clear mismatchInfoAll recept2clustAssignAll compTracksALTAll
            
            %The cut-off cluster size value we want is one less
            cutOffClust = cutOffClust - 1;
            
            %Write all remaining accumulated values into one file
            save(['valuesFor_',outFileID,'.mat'],...,
                'aggregStateAll','aggregStateFromR2CAll','clustHistoryMergedAll',...
                'eventTableAll','eventTable_monoAll','modelParamAll',...
                'rateAssocPerClustAll','rateDissocPerClustAll','rateAssocMonoAll',...
                'segmentStatAll','simParamAll','cutOffClust',...
                'elapsedTimeAll','assocStatsAll','collProbStatsAll',...
                '-v7.3');
            
            clear all
            
            fprintf('done.\n');
        end        
        
    end %if not input error
    
    
end %function
            
        