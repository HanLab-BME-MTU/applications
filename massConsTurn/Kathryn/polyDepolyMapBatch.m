function polyDepolyMapBatch
% POLYDEPOLYMAPBATCH: makes poly/depoly maps, time-averaged maps, and r/g colormaps of both

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-SET PARAMETERS - START
%
% firstWin/lastWin : Range of values to test for measurement window size
%                    (in pixels) in polyDepolyMap. Must be an odd number.
%                    The optimal value currently has to be determined
%                    empirically, but good starting values are 15/21.
%                    Please note: the window size is used to name the
%                    turnover directory, so in the case that you are NOT
%                    running poly/depoly, make these parameters correspond
%                    to those found in the directory names.
% doPolyDepolyCalc : 1 to run polyDepolyMap.m; 0 to skip.
% doNorm           : 1 to run image normalization within polyDepolyMap.m; 0
%                    to skip. Usually after you have run it at one window
%                    size, you won't need to repeat it. This parameter is
%                    irrelevant if doPolyDepolyCalc = 0.
%
% doTimeAvging     : 1 to run time averaging; 0 to skip.
% nFrms2Avg        : The number of frames to average. This parameter is
%                    irrelevant if doTimeAvging = 0.
% timeStepSize     : Frame interval for time averaging. This parameter is
%                    irrelevant if doTimeAvging = 0.
%
%
% makeRedGreenMaps : 1 to make red/green color maps of the data; 0 to skip.
%                    Data from raw and time-averaged .mat files are
%                    converted to color tif images.
% gamma            : Since most pixels will have a small value compared to
%                    a few with strong poly/depoly, you can boost the
%                    appearance of the signal by gamma correction.
%                    0<gamma<1 to boost low values.
% cMapLength       : Length of the color map.  Must be even. 128 should be
%                    fine in most cases.
% useGlobalMinMax  : 0 to normalize each movie's colormaps with its own values
%                    1 to normalize using values based on ALL the movies in the batch
%                    [userMin userMax] to use empirically derived values.

firstWin=17;  lastWin=17;
doPolyDepolyCalc=1;
doNorm=1;

doTimeAvging=1;
nFrms2Avg=5;
timeStepSize=1;

makeRedGreenMaps=1;
gamma=1;
cMapLength=128;
useGlobalMinMax=0; 

% USER-SET PARAMETERS - END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for winSize=firstWin:2:lastWin;

    topDir=uigetdir(pwd,'Please select top directory containing all the project directories');
    dirNames = dir(topDir);
    dirNames(1:2)=[]; % get rid of . and .. directories in list

    % get rid out of the list any element in dirNames that isn't a real
    % directory (e.g. .txt file, .mat file)
    countNonDirectories=zeros(length(dirNames)-2,1);
    for i=1:length(dirNames)
        countNonDirectories(i)=isdir([topDir filesep dirNames(i).name]);
    end
    dirNames(countNonDirectories==0)=[];
    nProj=length(dirNames);


    % run poly/depoly and time averaging if respective flags are 1
    batchMin=0; batchMax=0; % min and max values over the whole batch
    for i=1:nProj
        try
            imDir=[topDir filesep dirNames(i).name filesep 'images'];
            anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

            % run poly/depoly
            if doPolyDepolyCalc==1
                [runInfo]=polyDepolyMap(imDir,anDir,[],doNorm,winSize)
            else % or just get max/min info from previous run
                runInfo=load([anDir filesep 'turn_winSize_' num2str(winSize) filesep 'runInfo.mat']);
                runInfo=runInfo.runInfo;
            end

            % keep track of global min/max from all the movies
            batchMin=min(batchMin,runInfo.movieMin);
            batchMax=max(batchMax,runInfo.movieMax);

            % do time averaging
            if doTimeAvging==1
                polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize);
            end
        catch
            disp(['movie ' num2str(i) ' had an error during poly/depoly or time averaging:'])
            lasterr
        end

    end


    % now make red/green maps if flag is 1
    if makeRedGreenMaps==1
        if isequal(useGlobalMinMax,0)       % movie values
            minMax=[];
        elseif isequal(useGlobalMinMax,1)   % batch values
            minMax=[batchMin batchMax];
        else                                % user values
            minMax=useGlobalMinMax;
        end

        for i=1:nProj
            try
                cmDir=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'edge' filesep 'cell_mask'];
                turnDir=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'turn_winSize_' num2str(winSize)];
                mapDirSp=[turnDir filesep 'turnSpAvgDir'];
                mapDirTm=[turnDir  filesep 'turnTmAvgDir'];

                % get edge pixels to make border in red/green maps
                edgePix=load([turnDir filesep 'edgePix']);
                edgePix=edgePix.edgePix;

                % get r/g maps for spatial and temporal directories
                polyDepolyVisualizeActivity(cmDir,mapDirSp,edgePix,gamma,cMapLength,minMax);
                polyDepolyVisualizeActivity(cmDir,mapDirTm,edgePix,gamma,cMapLength,minMax);
            catch
                disp(['movie ' num2str(i) ' had an error during r/g map generation:'])
                lasterr
            end
        end
    end
    save([topDir filesep 'batchParams_winSize_' num2str(winSize)]);
end
