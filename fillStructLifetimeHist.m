function [data] = fillStructLifetimeHist(data);
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

% number of entries in exp
lens = length(data);

for i=1:lens
    
    % number of frames for this exp
    lenf = data(i).movieLength;
    
    % read lftInfo from the appropriate file if lftHist doesn't lareday
    % exist
    
    cvar = 1;
    if isfield(data,'lftHist')
        if ~isempty(data(i).lftHist)
            cvar = 0;
        end
    end
    
    if cvar==1

        path = data(i).source;
        lftpath = [path,'/LifetimeInfo'];
        od = cd;
        cd(lftpath);
        lftname = 'lftInfo.mat';
        loadfile = load(lftname);
        cd(od);
        lftInfo = loadfile.lftInfo;


        lftMat = lftInfo.Mat_lifetime;
        statMat =  lftInfo.Mat_status;

        [sx,sy] = size(lftMat);

        lftVec = nan*zeros(sx,1);

        % a trajectory is counted for the lifetime analysis if the status of 
        % the trajectory is ==1 and the value of the gaps is ==4
        for k=1:sx
            % current status vector
            cstat = nonzeros(statMat(k,:));
            % current lifetime vector
            clft = lftMat(k,:);
            if ( (min(cstat)==1) & (max(cstat)<5) )
                lftVec(k) = max(clft);
            end
        end    

        lftVec = lftVec(find(isfinite(lftVec)));
        data(i).lftHist = hist(lftVec,[1:lenf]);
        
    end % of if cvar==1
        

end % of for

end % of function
