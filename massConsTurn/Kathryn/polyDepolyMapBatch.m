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
% tmAvgAllFrmsAlso : 1 to average over the whole time series in addition
%                    to the run using nFrms2Avg. This parameter is
%                    irrelevant if doTimeAvging = 0.
%
% visualizeMovie   : 1 to make red/green color maps of the data; 0 to skip.
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
% mkMov            : 1 to make histogram and red/green map movies; 0 if not

winRange=[7:4:15];
doPolyDepolyCalc=0;
doNormInit=0;

doTimeAvging=0;
nFrms2Avg=5;
timeStepSize=1;
tmAvgAllFrmsAlso=0;

visualizeMovie=1;
gamma=1;
cMap='jet';
useMovValues=0; % 1 if you want to use the min/max values from the movie

% USER-SET PARAMETERS - END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

for polyLocFlag=1

    for winSize=winRange;
        
        if polyLocFlag==0 && winSize==winRange(1)
            doNorm=doNormInit;
        else
            doNorm=0;
        end
        
        for i=1:nProj

            imDir=[topDir filesep dirNames(i).name filesep 'images'];
            anDir=[topDir filesep dirNames(i).name filesep 'analysis-polyFlag' num2str(polyLocFlag)];

            try % run poly/depoly
                if doPolyDepolyCalc==1
                    disp(['Calculating polyDepoly maps: window size: ' num2str(winSize) ' pixels, Project ' num2str(i) ' of ' num2str(nProj)])
                    [polyDepoly,accumY,accumX,runInfo]=turnoverMap(imDir,anDir,[],doNorm,winSize,polyLocFlag)
                else % or just get max/min info from previous run
                    runInfo=load([anDir filesep 'turn_winL_' num2str(winSize) filesep 'runInfo.mat']);
                    runInfo=runInfo.runInfo;
                end
            catch
                disp(['error: Project ' num2str(i) ' during turnoverMap.m'])
                lasterr
            end


            try % make r/g maps from non-temporally averaged data
                if visualizeMovie==1
                    mapDir=[runInfo.turnDir filesep 'interp'];
                    polyDepolyVisualizeActivity(runInfo,mapDir,useMovValues,cMap,gamma,1);
                end
            catch
                disp(['error: Project ' num2str(i) ' during polyDepolyVisualizeActivity (1)'])
                lasterr
            end


            try % do time averaging
                if doTimeAvging==1
                    disp(['Calculating time-averaged maps: window size: ' num2str(winSize) ' pixels, Project ' num2str(i) ' of ' num2str(nProj)])
                    [runInfo]=polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize,[],[]);
                end
            catch
                disp(['error: Project ' num2str(i) ' during polyDepolyTimeAvg'])
                lasterr
            end


            try % make r/g maps from time-averaged data
                if visualizeMovie==1
                    mapDir=[runInfo.turnDir filesep 'timeAvg_' num2str(nFrms2Avg) '_' num2str(gamma) filesep 'mapMats'];
                    polyDepolyVisualizeActivity(runInfo,mapDir,useMovValues,cMap,gamma,1);
                end
            catch
                disp(['error: Project ' num2str(i) ' during polyDepolyVisualizeActivity (2)'])
                lasterr
            end
        end

    end
end
save([topDir filesep 'batchParams_winSize_' num2str(winSize)]);

