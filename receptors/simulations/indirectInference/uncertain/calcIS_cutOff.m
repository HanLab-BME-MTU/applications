function [isVals] = calcIS_cutOff(rootDir,numSims,outStr)
%CALCISCUTOFF calculates intermediate statistics for target and probe runs.
%This function is part of the indirect inference based model calibration
%framework.
%
%  INPUT:
%           rootDir:    location of outX directories with rates and also
%                       cluster statistics.
%           numSims:    number of simulations which corresponds to number
%                       of outX directories to process
%           outStr:     a string to be appended to the intermediate
%                       statistics output file name
%
%   OUTPUT:
%           isVals:     a structure with the following fields:
%                       1) theta - intermediate statistics vector
%                       2) intStatsMatrix - intermediate statistics matrix
%                       3) numDataPts - number of data points (sims) used
%                                       in calcuating theta
%                       4) largestClustSize - the largets cluster size
%                       5) cutOffClust - largest cluster size up to which
%                                        rates are calcuated
%                       6) mean_sd_rateAssocPerClust - mean and sd for
%                                                      association rates
%                       7) mean_sd_rateDissocPerClust - mean and sd for
%                                                       dissociation rates
%                       8) assocEventCount - number of association events
%                                            per cluster size
%                       9) dissocEventCount - number of dissociation events
%                                             per cluster size
%
%  Robel Yirdaw, 11/24/14
%

    %Intermediate statistic structure
    isVals = [];
    
    %Initialize read error flag and sim. index
    outputReadErr = 0;
    simIndx = 0;
    
    %Enable the following if using maximum cluster size restriction
    %MAX_CLUST_SIZE = 5;
    
    %If cluster statistics are located at a differet root, set the
    %following.
     sourceRoot = '';
    
    %Echo location and input parameters
    fprintf('\n======================================================================');
    fprintf('\ncalcIS at %s.',pwd);
    fprintf('\nrootDir is %s%s.',sourceRoot,rootDir);
    fprintf('\nnumSims is %d.',numSims);
    fprintf('\noutStr is %s.',outStr);
    fprintf('\n======================================================================');
    fprintf('\n');
        
    try 
        
        %Load clusterStats
        loadStruct = load([sourceRoot,rootDir,'/clusterStats_',outStr]);

        %The largest cluster size from target
        clusterStats = loadStruct.clusterStats;
        clear loadStruct
        largestClustSize = clusterStats.largestClustSize;
               
        %Will determine cutoff cluster size for association and
        %dissociation events.  Initialize to the maximum possible.
        cutOffClust = largestClustSize;

        %Will save calculated rates for individual runs
        rateAssocPerClustAll = NaN(numSims,largestClustSize);
        rateAssocMonoAll = NaN(numSims,largestClustSize);
        rateDissocPerClustAll = NaN(numSims,largestClustSize);
        %Will now collect event counts for each simulation
        dissocEventCount = NaN(numSims,largestClustSize);
        assocEventCount = NaN(numSims,largestClustSize);
        
        for simIndx=1:numSims
            fprintf('\nProcessing simulation #%d...',simIndx);
            
            %Construct path to specified location
            fullPathToRuns = [sourceRoot,rootDir,'/out',int2str(simIndx)];
            
            %Load output file containing rates
            loadStruct = load([fullPathToRuns,'/out',int2str(simIndx)]);
            
            %Cutoff cluster size for the current simulation. 10 or more
            %events must be present for association and dissociation events
            %for a given cluster size to be included
            currCutoffClust = find ( (loadStruct.eventTable(:,3) < 10) |...
                (loadStruct.eventTable(:,6) < 10),1,'first' ) - 1;
            
            rateAssocPerClustAll(simIndx,1:currCutoffClust) =...
                loadStruct.rateAssocPerClust(1:currCutoffClust);
            
            rateAssocMonoAll(simIndx,1:currCutoffClust) =...
                loadStruct.rateAssocMono(1:currCutoffClust);
            
            rateDissocPerClustAll(simIndx,1:currCutoffClust) =...
                loadStruct.rateDissocPerClust(1:currCutoffClust);
                        
            %Save event counts to determine the cutoff
            dissocEventCount(simIndx, 1:length(loadStruct.eventTable(:,3))) =...
                loadStruct.eventTable(:,3);
            assocEventCount(simIndx,1:length(loadStruct.eventTable(:,6))) =...
                loadStruct.eventTable(:,6);

            
            clear loadStruct currCutOffClust
            
            fprintf('done.');

        end %for each simulation output
        
        %Determine cutoff cluster size.
        %Only cluster sizes with at 10 or more association AND dissociation
        %events in at least 8 of the 10 simulations will be included. The
        %cutoff cluster is then the largest satisfying this criteria.

        cutOffClust = sum( sum(~isnan(rateAssocPerClustAll)) >= 8 );
        
    catch outputReadException
        fprintf('\n============================================================================');
        fprintf('\nException reading values (at simulation #%d).',simIndx);
        fprintf('\nRoot dir: %s \n\n',rootDir);       
        disp(outputReadException);
        fprintf('\n\n==========================================================================\n');        
        
        outputReadErr = 1;
    end
    
    if (~outputReadErr)
        
        try
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Need mean cluster density over iterations - numSim rows with
            %1:largestClustSize columns.
            
            clusterDensity_iterMean = transpose( reshape(mean(clusterStats.clusterDensity,2),...
                largestClustSize,numSims) );
            
            %Remove density values for erroneous sizes if necessary
            %clusterDensity_iterMean = clusterDensity_iterMean(:,1:min(largestClustSize,MAX_CLUST_SIZE));            

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Need normalized association rates.
            
            rateAssocMono = [];
            rateAssocClust = [];
            rateDissocClust = [];            
            mean_sd_rateAssocPerClust = [];
            mean_sd_rateDissocPerClust = [];
            
            %First normalize the rate of dimerization and save. There must
            %at lease be enough dimerization events
            if (cutOffClust > 0)
                rateAssocMono = rateAssocMonoAll(:,1)./(2*clusterDensity_iterMean(:,1));
                mean_sd_rateAssocPerClust(1:2,1:2) =...
                    [NaN NaN; nanmean(rateAssocMono) nanstd(rateAssocMono)];
            end
            
            if (cutOffClust > 1)
                rateAssocClust = rateAssocPerClustAll(:,2:cutOffClust)./...
                    repmat(clusterDensity_iterMean(:,1),1,(cutOffClust-1));
                mean_sd_rateAssocPerClust(3:cutOffClust+1,1:2) =...
                    [transpose(nanmean(rateAssocClust)) transpose(nanstd(rateAssocClust))];                
     
                rateDissocClust = rateDissocPerClustAll(:,2:cutOffClust);
                mean_sd_rateDissocPerClust =...
                    [NaN NaN; transpose(nanmean(rateDissocClust)) transpose(nanstd(rateDissocClust))];

            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Fill in intermediate statistics table and calculate the vector
            intStatsMatrix = [clusterDensity_iterMean rateAssocMono...
                rateAssocClust rateDissocClust];        

            %Finally calculate the theta vector
            theta = nanmean(intStatsMatrix);

            %Construct output structure
            isVals.theta = theta;
            isVals.intStatsMatrix = intStatsMatrix;
            isVals.numDataPts = sum(~isnan(intStatsMatrix));            
            isVals.largestClustSize = largestClustSize;
            isVals.cutOffClust = cutOffClust;
            isVals.mean_sd_rateAssocPerClust = mean_sd_rateAssocPerClust;
            isVals.mean_sd_rateDissocPerClust = mean_sd_rateDissocPerClust;
            isVals.assocEventCount = assocEventCount;
            isVals.dissocEventCount = dissocEventCount;

            %Output path and name
            outFile = [rootDir,'/isVals_',outStr,'_cutOff.mat'];
            
            %TEMP
            delete(outFile);
            %Write to file
            save(outFile,'isVals','-v7.3');
            
            fprintf('\n\nSaved intermediate statistics in %s.',outFile);
            fprintf('\n\nDone.');

        catch intStatsCreateException
            fprintf('\n============================================================================');
            fprintf('\nException creating intermediate statistics.\n');
            fprintf('\nRoot dir: %s \n\n',rootDir);
            disp(intStatsCreateException);
            fprintf('\n\n============================================================================\n');
            
        end
        
    end %if no read error

    
end %calcIS





