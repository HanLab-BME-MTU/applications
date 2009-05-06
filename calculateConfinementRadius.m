function [confRadius,trackCenter] = calculateConfinementRadius(tracks,probDim)
%calculateConfinementRadius outputs confinement radii for all tracks
%specified in input
%
%SYNOPSIS confinementRadius = calculateConfinementRadius(tracks,probDim);
%
%INPUT      tracks:            -- EITHER -- 
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%                           -- OR -- 
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%
%OUTPUT     confinementRadius:  Calculated the mean of the eigen vectors of 
%                               the variance-covariance matrix of this 
%                               track's positions 
%           trackCenter:        calculate the track's center
%
%USES       none
%
%created by DA Nunez April 2, 2009 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%estimate the confinement radius
for iTrack = 1:size(tracks,1)

    %get track coordinates
    xCoord = (tracks(iTrack,1:8:end))';
    yCoord = (tracks(iTrack,2:8:end))';
    zCoord = (tracks(iTrack,3:8:end))';
    xyzCoord = [xCoord yCoord zCoord];

    %find the eignevalues and eigenvectors of the variance-covariance
    %matrix of this track's positions
    eigenVal = eig(nancov(xyzCoord(:,1:probDim)));

    %calculate the track's confinement radius
    confRadius(iTrack,1) = sqrt( min(eigenVal) * (probDim + 2) );


    %calculate the track's center
    trackCenter(iTrack,:) = nanmean(xyzCoord(:,1:probDim));

end %for each track

end %of function