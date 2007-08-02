function dataStruct = makiGenerateTracks(dataStruct)
%MAKIGENERATETRACKS tracks kinetochores throughout a movie
%
%SYNOPSIS dataStruct = makiGenerateTracks(dataStruct)
%
%INPUT  dataStruct: dataStruct as in makiMakeDataStruct with at least the
%                   fields "dataProperties", "planeFit" & "initCoord".
%                   Optional. Loaded interactively if not input.
%
%OUTPUT dataStruct: Same as input, with added field "tracks"
%
%Khuloud Jaqaman, July 2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load dataStruct if not input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

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
    
%get number of features in each frame
if ~isfield(movieInfo,'num')
    for iTime = 1 : nTimepoints
        movieInfo(iTime).num = size(movieInfo(iTime).xCoord,1);
    end
end

%collect coordinates and their std in one matrix in each frame
if ~isfield(movieInfo,'allCoord')
    for iTime = 1 : nTimepoints
        movieInfo(iTime).allCoord = [movieInfo(iTime).xCoord ...
            movieInfo(iTime).yCoord movieInfo(iTime).zCoord];
    end
end

%calculate nearest neighbor distance for each feature in each frame
if ~isfield(movieInfo,'nnDist')

    for iTime = 1 : nTimepoints
        
        switch movieInfo(iTime).num

            case 0 %if there are no features

                %there are no nearest neighbor distances
                nnDist = zeros(0,1);

            case 1 %if there is only 1 feature

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000;

            otherwise %if there is more than 1 feature

                %compute distance matrix
                nnDist = createDistanceMatrix(movieInfo(iTime).allCoord(:,1:2:end),...
                    movieInfo(iTime).allCoord(:,1:2:end));

                %sort distance matrix and find nearest neighbor distance
                nnDist = sort(nnDist,2);
                nnDist = nnDist(:,3);

        end

        %store nearest neighbor distance
        movieInfo(iTime).nnDist = nnDist;

    end
    
end

%get tracking parameters
gapCloseParam = dataStruct.dataProperties.tracksParam.gapCloseParam;
costMatParam = dataStruct.dataProperties.tracksParam.costMatParam;
useLocalDensity = dataStruct.dataProperties.tracksParam.useLocalDensity;

%call tracker
try

    %track the kinetochores
    tracks = trackCloseGapsKalman(movieInfo,costMatParam,...
        gapCloseParam,[],useLocalDensity,0,3,0);

    %replace rotated coordinates in tracks by original ones if rotated 
    %coordinates were used for tracking
    if dataStruct.dataProperties.tracksParam.rotate == 1
        
        %go over all tracks
        for iTrack = 1 : length(tracks)
            
            %store rotated coordinates in another field
            tracks(iTrack).coordAmpRot = tracks(iTrack).tracksCoordAmpCG;
            
            %fetch the start and end time of this track
            startTime = tracks(iTrack).seqOfEvents(1,1);
            endTime = tracks(iTrack).seqOfEvents(2,1);
            
            %go over all frames where this track exists
            for iFrame =  startTime : endTime 
                
                %get the feature making up this track in this frame
                iFeature = tracks(iTrack).tracksFeatIndxCG(iFrame-startTime+1);
                
                %if there is a feature (not a gap)
                if iFeature ~= 0
                    
                    %replace coordiantes and their stds
                    tracks(iTrack).tracksCoordAmpCG(1,(iFrame-startTime)*8+1:...
                        (iFrame-startTime)*8+3) = dataStruct.initCoord(iFrame).allCoord(iFeature,1:3);
                    tracks(iTrack).tracksCoordAmpCG(1,(iFrame-startTime)*8+5:...
                        (iFrame-startTime)*8+7) = dataStruct.initCoord(iFrame).allCoord(iFeature,4:6);
                    
                end
                
            end
            
        end %(for iTrack = 1 : numTracks)
        
    end %(if dataStruct.dataProperties.tracksParam.rotate == 1)

catch
end

%store tracks in dataStruct
dataStruct.tracks = tracks;
