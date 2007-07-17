function dataStruct = makiGroupSisters(dataStruct)
%MAKIGROUPSISTERS groups sister kinetochores
%
% SYNOPSIS: dataStruct = makiGroupSisters(dataStruct)
%
% INPUT dataStruct: data structure as created by makiMakeDataStruct
%
% OUTPUT dataStruct: Same as input, but with field sisterList
%
% REMARKS The strategy is to first try and identify sister pairs among the
%           tracks that span at least 75% of the movie, and then to
%           successively try and group shorter tracks with the single
%           sisters.
%         Grouping is mainly based on the idea that within sister variance
%           is smaller than between sisters variance. Additionally, there
%           is a user-set distance cutoff above which pairs won't be linked
%         The code cannot handle merged/splitted tracks!
%
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 16-Jul-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defaults
goodTrackRatio = 0.75;
maxDist = 4; % max avg. distance in um, should be in dataProperties
robust = false; % whether to use robust statistics to estimate variance, mean


% test input
if nargin == 0 || isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

% check for tracks here in the future

%---

% read movieLength
nTimepoints = dataStruct.dataProperties.movieSize(4);
% read track statistics. This will work only if no merge/split, i.e. if
% there are only two events per track: a start and a finish
try
    trackStats = catStruct(3,'dataStruct.tracks.seqOfEvents');
catch
    error('makiGroupSisters cannot handle merging/splitting')
end


% check for relative length of tracks
relTrackLength = ...
    squeeze(trackStats(2,1,:)-trackStats(1,1,:)+1)./nTimepoints;
% select above 0.75
goodTracks = find(relTrackLength>=goodTrackRatio);
% count
nGoodTracks = length(goodTracks);

% preassign distance/variance matrix. Set zeros here for sparse later
[variances,distances,alignment] = deal(zeros(nGoodTracks));

% read normals to plane - takes time, thus outside the loop
        normals = catStruct(2,'dataStruct.planeFit.eigenVectors(:,1)')';

% loop through the good tracks, calculate for every pair median distance
% and variance
% DEBUG
figure
hold on
for jTrack = 1:nGoodTracks % loop cols

    % read index of track
    jIdx = goodTracks(jTrack);

    % read track coordinates
    colCoords = reshape(dataStruct.tracks(jIdx).tracksCoordAmpCG,8,[])';
    colCoords = colCoords(:,1:3);

    %     %DEBUG
    plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),'Color',extendedColors(jTrack))

    % read timepoints of the track
    colTime = (trackStats(1,1,jIdx):trackStats(2,1,jIdx))';
    % remove gaps
    colTime(all(isnan(colCoords),2)) = [];
    % remember indices into colCoords that correspond to the timepoints
    colIdx = colTime - trackStats(1,1,jIdx) + 1;

    for iTrack = jTrack+1:nGoodTracks % loop rows

        % read index of track
        iIdx = goodTracks(iTrack);

        % read track coordinates
        rowCoords = reshape(dataStruct.tracks(iIdx).tracksCoordAmpCG,8,[])';
        rowCoords = rowCoords(:,1:3);
        % read timepoints of the track
        rowTime = (trackStats(1,1,iIdx):trackStats(2,1,iIdx))';
        % remove gaps
        rowTime(all(isnan(rowCoords),2)) = [];
        % remember indices into colCoords that correspond to the timepoints
        rowIdx = rowTime - trackStats(1,1,iIdx) + 1;


        % find common time
        [commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

        % calculate distance (microns)
        distanceVector = colCoords(colIdx(ctColIdx),:) -...
            rowCoords(rowIdx(ctRowIdx),:);
        
        
        % calculate alignment:
        
        % retain commonTime-normals
        commonNormals = normals(commonTime,:);
        % cos(alpha) = dot(a,b)/(norm(a)*norm(b))
        [distance,distanceVectorN] = normList(distanceVector);
        alpha = acos(abs(dot(distanceVectorN',commonNormals')));
        % average alpha, rather than tan to be nice to pairs that will
        % align eventually
        meanAlpha = mean(alpha);
        
        
        
        
        % get average distance, variance
        % try robustStats
        if robust
            [rMean,rStd]=robustMean(distance);
        else
            rMean=mean(distance);
            rStd = std(distance);
        end
        distances(iTrack,jTrack) = rMean;
        variances(iTrack,jTrack) = rStd^2;
        distances(jTrack,iTrack) = rMean;
        variances(jTrack,iTrack) = rStd^2;
        % alignment: one up to about 30°, then increasing
        %[alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
        %    = deal(tan(meanAlpha)+1-meanAlpha);
        % start at one, increase from the beginning (2 at 45°)
        %[alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
        %    = deal(tan(meanAlpha)+1);
        % start at one, increase from the beginning (2 at 45°)
        %[alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
        %    = deal(tan(meanAlpha)+1);
        % start at one, increase from the beginning (2 at 30°)
        [alignment(jTrack,iTrack),alignment(iTrack,jTrack)]...
            = deal(2*sqrt(3)*tan(meanAlpha)+1);

    end
end

% cutoff distances

distCutoffIdx = distances>maxDist;
distances(distCutoffIdx) = 0;
variances(distCutoffIdx) = 0;
alignment(distCutoffIdx) = 0;

% make sparse cost matrix
%costMat = sparse(variances);
costMat = sparse(variances.*alignment);

% lap variances only
r2c = lap(costMat,-1,0,1);

% get cutoff for visualization

% find indices into matrix
ind=sub2ind([nGoodTracks nGoodTracks],1:nGoodTracks,double(r2c(1:nGoodTracks))');


% cutoff variances
[dummy,cutoff]=cutFirstHistMode(full(costMat(ind)),0);

% plot for one frame
t=[5,15,25,35];
plotGroupResults(t,r2c,nGoodTracks,goodTracks,dataStruct,costMat,cutoff)


% try: improve grouping by only considering 'good' variances - doesn't work
% well
%
% varCutoffIdx = variances>v*1.5;
% varOld = variances;
% variances(varCutoffIdx) = -1;
% r2c = lap(variances,-1,0,1);
%
% % plot for one frame
% t=5;
% plotGroupResults(t,r2c,nGoodTracks,goodTracks,dataStruct,variances)


% check for




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEBUG HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotGroupResults(tt,r2c,nGoodTracks,goodTracks,dataStruct,variances,cutoff)
% plot sister-links for one frame. Number indicates goodTrackIdx
% green: variance below cutoff
% blue : variance above cutoff
% red: track with no partner

if nargin < 7 || isempty(cutoff)
    cutoff = inf;
end
figure,
ntt = length(tt);
rows = ceil(sqrt(ntt));
cols = ceil(ntt/rows);
for p=1:ntt
    t = tt(p);
    [c2pg,c2pb,c2pr]=deal(nan(3*nGoodTracks,3));

    for i=nGoodTracks:-1:1
        idx1 = goodTracks(i);
        % check out of bounds
        if t<dataStruct.tracks(idx1).seqOfEvents(1,1) || t>dataStruct.tracks(idx1).seqOfEvents(2,1)
            c1 = nan(1,3);
        else
            % read via idx, not time
            tIdx = t-dataStruct.tracks(idx1).seqOfEvents(1,1)+1;
            c1=dataStruct.tracks(idx1).tracksCoordAmpCG((tIdx-1)*8+1:(tIdx-1)*8+3);
        end
        if r2c(i) > nGoodTracks
            c2=nan(1,3);
            v=NaN;
        else
            idx2 = goodTracks(r2c(i));
            if t<dataStruct.tracks(idx2).seqOfEvents(1,1) || t>dataStruct.tracks(idx2).seqOfEvents(2,1)
                c2 = nan(1,3);
            else
                % read via idx, not time
                tIdx = t-dataStruct.tracks(idx2).seqOfEvents(1,1)+1;
                c2=dataStruct.tracks(idx2).tracksCoordAmpCG((tIdx-1)*8+1:(tIdx-1)*8+3);
            end
            v = variances(i,r2c(i));
        end

        if v < cutoff
            c2pg(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        elseif isnan(v)
            c2pr(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        else
            c2pb(3*i-2:3*i,:) = [c1;c2;nan,nan,nan];
        end


    end
    subplot(rows,cols,p)
    plot3(c2pg(:,1),c2pg(:,2),c2pg(:,3),'*-g')
    hold on
    plot3(c2pb(:,1),c2pb(:,2),c2pb(:,3),'*-b')
    plot3(c2pr(:,1),c2pr(:,2),c2pr(:,3),'*-r')
    text(c2pg(1:3:end,1),c2pg(1:3:end,2),c2pg(1:3:end,3),num2str((1:nGoodTracks)'))
    text(c2pb(1:3:end,1),c2pb(1:3:end,2),c2pb(1:3:end,3),num2str((1:nGoodTracks)'))
    text(c2pr(1:3:end,1),c2pr(1:3:end,2),c2pr(1:3:end,3),num2str((1:nGoodTracks)'))
    grid on
end


function distance = plotTrackDistance(iIdxT,jIdxT,goodTracks,dataStruct,trackStats)
% plots distance between two tracks vs. time; two connected tracks

jIdx = goodTracks(jIdxT);
iIdx = goodTracks(iIdxT);

[rowCoords,rowTime,rowIdx] = trackData(iIdx,dataStruct,trackStats);
[colCoords,colTime,colIdx] = trackData(jIdx,dataStruct,trackStats);
% find common time
[commonTime,ctColIdx,ctRowIdx] = intersect(colTime,rowTime);

% calculate distance (microns)
distance = sqrt(sum((colCoords(colIdx(ctColIdx),:) -...
    rowCoords(rowIdx(ctRowIdx),:)).^2,2));

figure('Name',sprintf('tracks %i (r) & %i (b)',iIdxT,jIdxT))

subplot(1,2,1)
plot(commonTime,distance)
ylim([0,max(5,max(distance))])

% also plot 3d track
subplot(1,2,2)
plot3(colCoords(:,1),colCoords(:,2),colCoords(:,3),'b')
hold on
plot3(rowCoords(:,1),rowCoords(:,2),rowCoords(:,3),'r')
cc = reshape([colCoords(colIdx(ctColIdx),:),...
    rowCoords(rowIdx(ctRowIdx),:),nan(length(distance),3)]',3,[])';
plot3(cc(:,1),cc(:,2),cc(:,3),'g')
plot3(cc(1:2,1),cc(1:2,2),cc(1:2,3),'*')

function [coords,time,coordIdx] = trackData(idx,dataStruct,trackStats)

% read track coordinates for j
coords = reshape(dataStruct.tracks(idx).tracksCoordAmpCG,8,[])';
coords = coords(:,1:3);
% read timepoints of the track
time = (trackStats(1,1,idx):trackStats(2,1,idx))';
% remove gaps
time(all(isnan(coords),2)) = [];
% remember indices into colCoords that correspond to the timepoints
coordIdx = time - trackStats(1,1,idx) + 1;






