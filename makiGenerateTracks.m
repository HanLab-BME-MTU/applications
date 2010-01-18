function dataStruct = makiGenerateTracks(dataStruct)
%MAKIGENERATETRACKS tracks kinetochores throughout a movie
%
%SYNOPSIS dataStruct = makiGenerateTracks(dataStruct)
%
%INPUT  dataStruct: dataStruct as in makiMakeDataStruct with the
%                   fields "dataProperties", "initCoord" & "planeFit".
%                   Field "planeFit" can be empty.
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

%if rotated coordinates are to be used ...
if dataStruct.dataProperties.tracksParam.rotate == 1 && ~isempty(dataStruct.planeFit)
    
    %get the rotated coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.planeFit(iTime).rotatedCoord;
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
    end
    
else %if the original coordinates are to be used
    
    %get the original coordinates
    for iTime = 1 : nTimepoints
        allCoord = dataStruct.initCoord(iTime).allCoord;
        movieInfo(iTime).xCoord = [allCoord(:,1) allCoord(:,4)];
        movieInfo(iTime).yCoord = [allCoord(:,2) allCoord(:,5)];
        movieInfo(iTime).zCoord = [allCoord(:,3) allCoord(:,6)];
        movieInfo(iTime).amp = dataStruct.initCoord(iTime).amp;
    end
    
    %if user requested compensating for COM translation
    if dataStruct.dataProperties.tracksParam.rotate == 0
        
        %calculate the center of mass in each frame
        centerOfMass = zeros(nTimepoints,3);
        for iTime = 1 : nTimepoints
            centerOfMass(iTime,:) = [mean(movieInfo(iTime).xCoord(:,1)) ...
                mean(movieInfo(iTime).yCoord(:,1)) mean(movieInfo(iTime).zCoord(:,1))];
        end
        
        %shift coordinates by center of mass to make the origin in each frame
        %at its center of mass
        for iTime = 1 : nTimepoints
            movieInfo(iTime).xCoord(:,1) = movieInfo(iTime).xCoord(:,1) - centerOfMass(iTime,1);
            movieInfo(iTime).yCoord(:,1) = movieInfo(iTime).yCoord(:,1) - centerOfMass(iTime,2);
            movieInfo(iTime).zCoord(:,1) = movieInfo(iTime).zCoord(:,1) - centerOfMass(iTime,3);
        end
        
    end
    
end

%if any coordinate is exactly 0, add to it eps in order to avoid getting
%into a mess with sparse format
for iTime = 1 : nTimepoints
    coordTmp = movieInfo(iTime).xCoord;
    movieInfo(iTime).xCoord(coordTmp==0) = eps;
    coordTmp = movieInfo(iTime).yCoord;
    movieInfo(iTime).yCoord(coordTmp==0) = eps;
    coordTmp = movieInfo(iTime).zCoord;
    movieInfo(iTime).zCoord(coordTmp==0) = eps;
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

            case {1,2} %if there are only 1 or 2 features

                %assign nearest neighbor distance as 1000 pixels (a very big
                %number)
                nnDist = 1000*ones(movieInfo(iTime).num,1);
                
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

%get kinetochore classification in each frame if available
if ~isempty(dataStruct.planeFit) %if the plane fit has been done
    
    %extract planeFit field from structure
    planeFit = dataStruct.planeFit;
    
    %assign kinetochore types per frame
    for iTime = 1 : nTimepoints
        kinType = zeros(movieInfo(iTime).num,1); %inlier
        kinType(planeFit(iTime).unalignedIdx) = 1; %unaligned
        kinType(planeFit(iTime).laggingIdx) = 2; %lagging
        movieInfo(iTime).kinType = kinType;
    end
    
else %if not
    
    %treat all kinetochores as inliers
    for iTime = 1 : nTimepoints
        movieInfo(iTime).kinType = zeros(movieInfo(iTime).num,1);
    end

end

%get tracking parameters
gapCloseParam = dataStruct.dataProperties.tracksParam.gapCloseParam;
costMatrices = dataStruct.dataProperties.tracksParam.costMatrices;
kalmanFunctions = dataStruct.dataProperties.tracksParam.kalmanFunctions;

%call tracker
try
    
    %track the kinetochores
    tracks = trackCloseGapsKalmanSparse(movieInfo,...
        costMatrices,gapCloseParam,kalmanFunctions,3,0,0);

    %replace the coordinate used for tracking (whether the rotated
    %coordinates or the original coordinates shifted by center of mass) by
    %the original coordinates
    for iTrack = 1 : length(tracks)

        %store coordinates used for tracking in another field
        tracks(iTrack).coordAmp4Tracking = tracks(iTrack).tracksCoordAmpCG;

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

catch
    disp('error in tracking')
end

%store tracks in dataStruct
dataStruct.tracks = tracks;
