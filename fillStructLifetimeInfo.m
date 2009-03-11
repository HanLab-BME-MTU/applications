function [exp]=fillStructLifetimeInfo(exp);
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

% number of entries in exp
lens = length(exp);

for i=1:lens
    % number of frames for this exp
    lenf = exp(i).movieLength;
    
    % first determine if lifetime data already exists
    lftVar = 0;
    
    if ~isfield(exp,'lftInfo')
        lftVar = 1;
    else
        if isempty(exp(i).lftInfo)
            lftVar = 1;
        end
    end
    
    % if lifetime data still needs to be determined
    if lftVar == 1
        % read trackInfo from exp if it exists as a field there; if not,
        % read it from the appropriate file
        readVar = 0;
        if isfield(exp,'trackInfo')
            readfield = exp(i).trackInfo;
            if ~isempty(readfield)
                trackInfo = readfield;
            else
                readVar = 1;
            end
        else
            readVar = 1;
        end
        
        if readVar == 1
            path = exp(i).source;
            trackpath = [path,'/TrackInfoMatrices'];
            od = cd;
            cd(trackpath);
            % current trackInfo name depends on whether it's a resampled
            % movie
            trackname = 'trackInfo.mat';
            
            loadfile = load(trackname);
            cd(od);
            
            trackInfo = loadfile.trackInfo;
        end % of if trackInfo has to be loaded
            
        
        
        [lftMat,statMat,xMat,yMat,disappMat,lftVec]=findLifetimesStatusSimple(trackInfo);
        exp(i).lftInfo.Mat_lifetime = lftMat;
        exp(i).lftInfo.Mat_status = statMat;
        exp(i).lftInfo.Mat_xcoord = xMat;
        exp(i).lftInfo.Mat_ycoord = yMat;
        exp(i).lftInfo.Mat_disapp = disappMat;
        
        exp(i).lftVec = lftVec;
        exp(i).lftHist = hist(lftVec,[1:lenf]);
        
%     end % of for

end % of for

end % of function
