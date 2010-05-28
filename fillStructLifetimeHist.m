function [data] = fillStructLifetimeHist(data, overwrite)
% fill experiment structure with lifetime histogram based on the lftInfo
% file saved at the specified directory
% FILLSTRUCTLIFETIMEINFO calculates lifetime from trackInfo matrix and
% fills them into exp structure
%
% SYNOPSIS [exp]=fillStructLifetimeInfo(exp);
%
% INPUT     data    :    experiment structure, which has to contain the
%                       fields
%                       .source
%                       .movieLength
%                       source is the path to the data location; at this
%                       location, the function reads the lftInfo
%                       from a folder called LifetimeInfo
% OUTPUT    exp         creates a new field
%                       .lftHist   
% REMARKS 
%
% Dinah Loerke, last modified Mar 2008
% Francois Aguet, 01/21/2010

if nargin<2
    overwrite = 0;
end

for i = 1:length(data)
    
    if ~isfield(data,'lftHist') || isempty(data(i).lftHist) || overwrite

        load([data(i).source 'LifetimeInfo' filesep 'lftInfo.mat']);

        lftMat = full(lftInfo.Mat_lifetime);
        statMat = full(lftInfo.Mat_status);

        sx = size(lftMat,1);
        lftVec = NaN(sx,1);

        % a trajectory is counted for the lifetime analysis if the status of 
        % the trajectory is ==1 and the value of the gaps is ==4
        for k = 1:sx
            % current status vector
            cstat = nonzeros(statMat(k,:));
            if ( (min(cstat)==1) && (max(cstat)<5) )
                lftVec(k) = max(lftMat(k,:));
            end
        end    

        lftVec = lftVec(isfinite(lftVec));
        data(i).lftHist = hist(lftVec, 1:data(i).movieLength);
    end
end