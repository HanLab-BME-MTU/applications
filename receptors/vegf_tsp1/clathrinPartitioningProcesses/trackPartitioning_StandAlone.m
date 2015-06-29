function [partitionResult] = trackPartitioning_StandAlone(tracks, mask, xMax, yMax, isSingleFrame)
%Determines the paritioning fraction of tracks into the mask
%
%SYNOPSIS [partitionResult] = trackPartitioning_StandAlone(tracks, mask, xMax, yMax, isSingleFrame)
%
%INPUT
%   tracks          : contains information about the tracks to be analyzed
%   mask            : logical array used to analyze how tracks partition
%                     into 'true' elements of the mask.
%   xMax            : size of the images (in image coordinates)
%   yMax            : size of the images
%   isSingleFrame   : logical is number of frames in mask == 1
%
%OUTPUT
%   partitionResult : contains partitioing information of the tracks
%       .nFramesTot     : the total number of frames the track was seen in
%       .nFramesIn      : the total number of frames the track was within
%                         the mask
%       .partitionFrac  : the partitioning fraction. 0 = the track was
%                         never found within mask. 1 = the track was always
%                         found within the mask.
%
%Notes
%   Mask array is in plot coordinate system [y, x]
%
%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
%standard process input
ip.addRequired('tracks', @isstruct);
ip.addRequired('mask', @islogical);
ip.addRequired('xMax', @isnumeric);
ip.addRequired('yMax', @isnumeric);
ip.addRequired('isSingleFrame', @islogical);
ip.parse(tracks, mask, xMax, yMax, isSingleFrame);
%% Partitioning Analysis
%calls function that does partititoning analysis
partitionResult = arrayfun(@partition, tracks);
%% Nested function
    function [result] = partition(data)
        [nSubtracks, nCoord] = size(data.tracksCoordAmpCG);
        for iSubtracks = 1:nSubtracks
            nFramesTot = 0;
            nFramesIn = 0;
            frame = data.seqOfEvents(1,1);
            for iCoord = 1:8:nCoord
                x = round(data.tracksCoordAmpCG(iSubtracks, iCoord));
                y = round(data.tracksCoordAmpCG(iSubtracks, iCoord+1));
                if x>0 && x<=xMax && y>0 && y<=yMax
                    nFramesTot = nFramesTot + 1;
                    if isSingleFrame
                        if mask(y,x)
                            nFramesIn = nFramesIn +1;
                        end
                    else
                        if mask(y,x,frame)
                            nFramesIn = nFramesIn +1;
                        end
                    end
                end
                frame = frame +1;
            end
            result.nFramesTot(iSubtracks) = nFramesTot;
            result.nFramesIn(iSubtracks) = nFramesIn;
            result.partitionFrac(iSubtracks) = nFramesIn / nFramesTot;
        end
    end


end

