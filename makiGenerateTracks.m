function dataStruct = makiGenerateTracks(dataStruct)
%MAKIGENERATETRACKS tracks kinetochores throughout a movie
%
%SYNOPSIS dataStruct = makiGenerateTracks(dataStruct)
%
%INPUT dataStruct: dataStruct as in makimakeDataStruct with at least the
%                   fields
%                     .dataProperties
%                     .planeFit
%                     .initCoord
%
%OUTPUT dataStruct.tracks
%
%Khuloud Jaqaman, July 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of time points in movie
nTimepoints = dataStruct.dataProperties.movieSize(4);

%get kinetochore coordinates and amplitude
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]),...
    nTimepoints,1);
if dataStruct.dataProperties.tracksParam.rotate == 1 %rotated coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.planeFit(iTime).rotatedCoord;
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
    end
else %original coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.initCoord(iTime).allCoord;
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
    end
end
    
%get tracking parameters
gapCloseParam = dataStruct.dataProperties.tracksParam.gapCloseParam;
costMatParam = dataStruct.dataProperties.tracksParam.costMatParam;
useLocalDensity = dataStruct.dataProperties.tracksParam.useLocalDensity;

%call tracker
try
    tracks = trackCloseGapsKalman(movieInfo,costMatParam,...
        gapCloseParam,[],useLocalDensity,0,3,0);
catch
end

%store tracks in dataStruct
dataStruct.tracks = tracks;
