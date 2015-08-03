function calculateSimSetQuants(rootDir,numSims,outStr)
%CALCULATESIMSETQUANTS calculates clusterStats by calling
%calcClusterStatsFromCompTracks. This function is part of the indirect 
%inference based model calibration framework.
%
%   INPUT:
%           rootDir:    location to outX directories with compTracksALT
%           numSims:    number of simulations which corresponds to number
%                       of outX directories/compTracksALT to process
%           outStr:     a string to be appended to the cluster statistics
%                       output file name
%
%   OUTPUT:
%           cluster statistics are saved to a mat file called
%           clusterStats_outStr.
%           
%   Robel Yirdaw, November 2014
%

    %Initialize file read error flag
    outputReadErr = 0;
    %Echo input parameters
    fprintf('\n======================================================================');
    fprintf('\nAt %s.',pwd);
    fprintf('\nrootDir is %s.',rootDir);
    fprintf('\nnumSims is %d.',numSims);
    fprintf('\noutStr is %s.',outStr);
    fprintf('\n======================================================================');
    fprintf('\n');
    
    try        
        %Will need compTracksALT
        compTracksALTAll = cell(numSims,1);
                
        for simIndx=1:numSims
            fprintf('\nProcessing simulation #%d...',simIndx);
            
            %Construct path to current specified location
            fullPathToRuns = [rootDir,'/out',int2str(simIndx)];

            %Load compTracksALT
            loadStruct = load([fullPathToRuns,'/compTracksALT',int2str(simIndx)]);
            compTracksALTAll{simIndx} = loadStruct.compTracksALT;
            
            clear loadStruct
            
            fprintf('done.');

        end %for each simulation output
        
    catch outputReadException
        fprintf('\n============================================================================');
        fprintf('\nException reading values (at simulation #%d).',simIndx);
        fprintf('\nRoot dir: %s \n\n',rootDir);       
        disp(outputReadException);
        fprintf('\n\n============================================================================\n');        
        
        outputReadErr = 1;
    end
    
    if (~outputReadErr)

        fprintf('\n\nCalculate clusterStats.\n');        
        
        try

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %compTracksALT for repeat runs loaded.
            %Next determine cluster statistics from compTracks.
            observeSideLen = 25;
            probDim = 2;
            clusterStats = calcClusterStatsFromCompTracks(compTracksALTAll,observeSideLen,probDim);
            %Save clusterStats
            save([rootDir,'/clusterStats_',outStr,'.mat'],'clusterStats','-v7.3');
            
            fprintf('\nSaved clusterStats.');

            fprintf('\n\nDone.');
            
            clear clusterStats
            
        catch calcException
            fprintf('\n============================================================================');
            fprintf('\nException performing calculations.\n');
            fprintf('\nRoot dir: %s \n\n',rootDir);
            disp(calcException);
            fprintf('\n\n============================================================================\n');
            
        end
        
    end %if no read error

    
end %calcSSQ





