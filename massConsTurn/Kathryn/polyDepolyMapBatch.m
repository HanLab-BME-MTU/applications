function polyDepolyMapBatch
% POLYDEPOLYMAPBATCH: makes poly/depoly maps, time-averaged maps, and r/g colormaps of both

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER-SET PARAMETERS - START


doNorm=0; % run image normalization

chooseROI=0; % user picks ROI for each project (do after normalization)

doVecOutlier=0; % get vector coverage masks by flow outlier detection
thresh=50;
nNeighborLowerLimit=1;
dirDev=45;

doTurnover=1; % run turnoverMap
winRange=[7]; %:4:15];
methodStr='frame0grid'; % 'frame0grid', 'frame1grid', or 'CAW'

doTimeAvging=1; % time averaging
nFrms2Avg=4;
timeStepSize=1;

doVisualize=1; % make tifs for turnoverMap and time-averaged data
gamma=1;
cMap='iso'; % 'jet' or 'iso'
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

if doNorm==1
    for i=1:nProj
        runInfo.imDir=[topDir filesep dirNames(i).name filesep 'images'];
        runInfo.anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

        normImgSeries(runInfo,0);
    end
end

if chooseROI==1
    for i=1:nProj
        runInfo.imDir=[topDir filesep dirNames(i).name filesep 'images'];
        runInfo.anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

        polyDepolyChooseRoi(runInfo);
    end
end

if doVecOutlier==1
    for i=1:nProj
        runInfo.imDir=[topDir filesep dirNames(i).name filesep 'images'];
        runInfo.anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

        runInfo.thresh=thresh;
        runInfo.nNeighborLowerLimit=nNeighborLowerLimit;
        runInfo.dirDev=dirDev;

        removeVectorOutliers(runInfo);
    end
end



for winL=winRange;
    for i=1:nProj

        imDir=[topDir filesep dirNames(i).name filesep 'images'];
        anDir=[topDir filesep dirNames(i).name filesep 'analysis'];

        temp=load([anDir filesep 'runInfo.mat']);
        runInfo=temp.runInfo;

        try % run poly/depoly
            if doTurnover==1
                disp(['Calculating polyDepoly maps: window size: ' num2str(winL) ' pixels, Project ' num2str(i) ' of ' num2str(nProj)])
                if isequal(methodStr,'frame0grid') || isequal(methodStr,'frame1grid')
                    turnoverMap(runInfo,[],winL,methodStr);
                elseif isequal(methodStr,'CAW')
                    turnoverMap_CAW(runInfo,[],winL);
                end
            end
        catch
            disp(['error: Project ' num2str(i) ' during turnoverMap.m'])
            lasterr
        end
        
        try % get runInfo from specific turnover directory
            turnDir=[anDir filesep 'turn_' methodStr filesep 'turn_winL_' num2str(winL)];
            temp=load([turnDir filesep 'runInfo.mat']);
            runInfo=temp.runInfo;
        catch
            disp('error: cannot do time-averaging and/or visualization without running turnoverMap first')
            lasterr
        end
        
        try % make r/g maps from non-temporally averaged data
            if doVisualize==1
                mapDir=[runInfo.turnDir filesep 'interp'];
                polyDepolyVisualizeActivity(runInfo,mapDir,useMovValues,cMap,gamma,0);
            end
        catch
            disp(['error: Project ' num2str(i) ' during polyDepolyVisualizeActivity (1)'])
            lasterr
        end


        try % do time averaging
            if doTimeAvging==1
                disp(['Calculating time-averaged maps: window size: ' num2str(winL) ' pixels, Project ' num2str(i) ' of ' num2str(nProj)])
                polyDepolyTimeAvg(runInfo,nFrms2Avg,timeStepSize,[],[]);
            end
        catch
            disp(['error: Project ' num2str(i) ' during polyDepolyTimeAvg'])
            lasterr
        end


        try % make r/g maps from time-averaged data
            if doVisualize==1
                mapDir=[runInfo.turnDir filesep 'timeAvg_' num2str(nFrms2Avg) '_' num2str(gamma) filesep 'mapMats'];
                polyDepolyVisualizeActivity(runInfo,mapDir,useMovValues,cMap,gamma,0);
            end
        catch
            disp(['error: Project ' num2str(i) ' during polyDepolyVisualizeActivity (2)'])
            lasterr
        end
    end

end

