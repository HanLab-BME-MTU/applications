function polyDepolyMapBatch
% makes poly/depoly maps, does time averaging, and makes r/g colormaps of
% both raw and time-averaged maps

%%%%%%%%%%%%%%%%%%%%%
% user can change these
doNorm=1; % 1 if you want to do image normalization for sure; if 0 it only does it if it hasn't been run before
gamma=1; % should be <1 if you want to boost low values
cMapLength=128; % colormap length; 128 should be fine. must be even.
nFrms2Avg=5; % for time averaging, the number of frames to average
timeStepSize=1; % frame interval for time averaging
useGlobalMaxMin=0; % 1 if you want to compare r/g maps between movies; 0 if you want each movie to have max/min red/green
doPolyDepolyCalc=1; % 1 if you need to run polyDepolyMap first; 0 if you just want to play with r/g maps
winSize=15;
%%%%%%%%%%%%%%%%%%%%%

topDir=uigetdir(pwd,'Please select top directory containing all the project directories');
dirNames = dir(topDir);

% here we remove entries for . and .. and the txt file at the end of the
% control cell list
dirNames(1:2)=[];
dirNames(end)=[];

nProj=length(dirNames);

% these will keep track of min and max values over the whole batch
batchMin=0; batchMax=0;
if doPolyDepolyCalc==1
    for i=1:nProj
        imDir=[topDir filesep dirNames(i).name filesep 'images'];
        anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

        % run poly/depoly
        [runInfo]=polyDepolyMap(imDir,anDir,[],doNorm,winSize)

        % do time averaging
        polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize)

        % keep track of global min/max from all the movies
        batchMin=min(batchMin,runInfo.movieMin);
        batchMax=max(batchMax,runInfo.movieMax);
    end
end

for i=1:nProj % get r/g maps all normalized to same values %%% analysis\turn\turnSpAvgDir
    cmDir=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'edge' filesep 'cell_mask'];
    mapDirSp=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'turn' filesep 'turnSpAvgDir'];
    mapDirTm=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'turn' filesep 'turnTmAvgDir'];

    % get edge pixels to make border in red/green maps
    turnDir=[topDir filesep dirNames(i).name filesep 'analysis' filesep 'turn'];
    edgePix=load([turnDir filesep 'edgePix'])
    edgePix=edgePix.edgePix;

    if useGlobalMaxMin==1 % use same values for all movies to make r/g maps (so can compare visually)
        polyDepolyVisualizeActivity(cmDir,mapDirSp,edgePix,gamma,cMapLength,[batchMin batchMax]); % each individual frame
        polyDepolyVisualizeActivity(cmDir,mapDirTm,edgePix,gamma,cMapLength,[batchMin batchMax]); % each time-averaged frame (same normalization)
    else % use values appropriate to each movie to make r/g maps (so each looks good, but you can't compare movies)
        polyDepolyVisualizeActivity(cmDir,mapDirSp,edgePix,gamma,cMapLength); % each individual frame
        polyDepolyVisualizeActivity(cmDir,mapDirTm,edgePix,gamma,cMapLength); % each time-averaged frame (same normalization)
    end
end

save([topDir filesep 'batchMinMax'],'batchMin','batchMax','gamma','cMapLength','dirNames','nFrms2Avg','timeStepSize');
