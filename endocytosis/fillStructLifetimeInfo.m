function [exp] = fillStructLifetimeInfo(exp)
% fill experiment structure with lifetime info based on the existing
% trackInfo
% FILLSTRUCTLIFETIMEINFO calculates lifetime from trackInfo matrix and
% fills them into exp structure
%
% SYNOPSIS [exp]=fillStructLifetimeInfo(exp);
%
% INPUT     exp    :    experiment structure, which has to contain the
%                       fields
%                       .source
%                       .movieLength
%                       source is the path to the data location; at this
%                       location, the function reads the trackInfo
%                       from a folder called TrackInfoMatrices
% OUTPUT    exp         creates a number of new fields, which are
%                       .lftInfo.Mat_lifetime = lftMat;
%                       .lftInfo.Mat_status = statMat;
%                       .lftInfo.Mat_xcoord = xMat;
%                       .lftInfo.Mat_ycoord = yMat;
%                       .lftInfo.Mat_disapp = disappMat;
%                       .lftVec = lftVec;
%                       .lftHist = hist(lftVec,[1:lenf]);
% REMARKS
%
% Dinah Loerke, last modified Jan 2008
% Francois Aguet, December 2009

for i = 1:length(exp)  
    
    if ~isfield(exp(i), 'lftInfo') || isempty(exp(i).lftInfo)
        if ~isfield(exp(i), 'trackInfo') || isempty(exp(i).trackInfo)
            load([exp(i).source 'TrackInfoMatrices' filesep 'trackInfo.mat']);
            exp(i).trackInfo = trackInfo;
        end
        [lftMat, statMat, xMat, yMat, disappMat, lftVec] = findLifetimesStatusSimple(exp(i).trackInfo);
        exp(i).lftInfo.Mat_lifetime = lftMat;
        exp(i).lftInfo.Mat_status = statMat;
        exp(i).lftInfo.Mat_xcoord = xMat;
        exp(i).lftInfo.Mat_ycoord = yMat;
        exp(i).lftInfo.Mat_disapp = disappMat;
        exp(i).lftVec = lftVec;
        exp(i).lftHist = hist(lftVec, 1:exp(i).movieLength);
    end    
end