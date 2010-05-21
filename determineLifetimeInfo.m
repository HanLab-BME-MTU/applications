function [data] = determineLifetimeInfo(data, overwrite)
% determine lifetime info data for all entries in the data structure
%
% SYNOPSIS  [data]=determineLifetimeInfo(data);
%
% INPUT     data    :    experiment structure, which has to contain a field
%                       .source, which is the path to the data location; at
%                       this location, the function reads the trackInfo
%                       from a folder called TrackInfoMatrices
% OUTPUT    data
% REMARKS
%
% Dinah Loerke, last modified Feb 2008
% Francois Aguet, December 2009

if (nargin < 2)
    overwrite = 0;
end

for i=1:length(data)
    lftPath = [data(i).source 'LifetimeInfo' filesep 'lftInfo.mat'];    
    if ~(exist(lftPath, 'file')==2) || (overwrite==1)
        fprintf('Generating lifetimes for movie %d of %d\n', i, length(data));
        
        trackedFeaturesPath = [data(i).source 'TrackInfoMatrices' filesep 'trackedFeatures.mat'];
        tracksFinalPath = [data(i).source 'TrackInfoMatrices' filesep 'tracksFinal.mat'];
        
        if (exist(trackedFeaturesPath, 'file')==2)
            tfile = load(trackedFeaturesPath);
            trackInfo = tfile.trackedFeatureInfo;
        elseif (exist(tracksFinalPath, 'file')==2)
            tfile = load(tracksFinalPath);
            if isfield(tfile, 'tracksFinal')
                trackInfo = convStruct2MatNoMS(tfile.tracksFinal);
            else
                trackInfo = [];
            end
        else
            trackInfo = [];
        end
        
        if ~isempty(trackInfo)
            [lftMat, statMat, xMat, yMat, disappMat] = determineLifetimeStatus(trackInfo);
            lftInfo.Mat_xcoord = xMat;
            lftInfo.Mat_ycoord = yMat;
            lftInfo.Mat_status = statMat;
            lftInfo.Mat_lifetime = lftMat;
            lftInfo.Mat_disapp = disappMat;
            
            lftDir = [data(i).source 'LifetimeInfo'];
            if ~(exist(lftDir, 'dir')==7)
                mkdir(lftDir);
            end
            save([lftDir filesep 'lftInfo.mat'], 'lftInfo');
        end
    end
end