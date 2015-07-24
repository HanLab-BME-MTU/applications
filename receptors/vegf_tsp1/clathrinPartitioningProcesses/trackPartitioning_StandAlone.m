function [partitionResult] = trackPartitioning_StandAlone(tracks, mask, ROIMask, xMax, yMax, isSingleFrame, varargin)
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
ip.KeepUnmatched = true;
%standard process input
ip.addRequired('tracks', @isstruct);
ip.addRequired('mask', @islogical);
ip.addRequired('xMax', @isnumeric);
ip.addRequired('yMax', @isnumeric);
ip.addRequired('isSingleFrame', @islogical);
ip.addParameter('scrambleTracks', false, @(x) islogical(x) || isnumeric(x));
ip.parse(tracks, mask, xMax, yMax, isSingleFrame, varargin{:});
scrambleTracks = ip.Results.scrambleTracks;
%% Partitioning Analysis
%calls function that does partititoning analysis
partitionResult = arrayfun(@partition, tracks);
%% Nested function
    function [result] = partition(data)
        [nSubtracks, nCoord] = size(data.tracksCoordAmpCG);
        for iSubtracks = 1:nSubtracks
            %determine mean loaction-----------------------
            nLoc = nCoord / 8;
            x_all = zeros(1, nLoc);
            y_all = zeros(1, nLoc);
            indx = 1;
            for iCoord = 1:8:nCoord
                x_all(indx) = data.tracksCoordAmpCG(iSubtracks, iCoord);
                y_all(indx) = data.tracksCoordAmpCG(iSubtracks, iCoord+1);
                indx = indx + 1;
            end
            if scrambleTracks
                x_mean = mean(x_all(~isnan(x_all)));
                y_mean = mean(y_all(~isnan(y_all)));
                x_all = x_all - x_mean + rand() * xMax + 0.5;
                y_all = y_all - y_mean + rand() * yMax + 0.5;
            end
            %determine partition fraction-------------------
            nFramesTot = 0;
            nFramesIn = 0;
            frame = data.seqOfEvents(1,1);
            for iLoc = 1:nLoc
                x = round(x_all(iLoc));
                y = round(y_all(iLoc));
                if x>0 && x<=xMax && y>0 && y<=yMax && ROIMask(y, x)
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

