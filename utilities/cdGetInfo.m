function idlistList = cdGetInfo(varargin)
%CDGETINFO allows extracting specific information from idlists
%
% SYNOPSIS: out = cdGetInfo(identifier1, identifier2,...)
%
% INPUT identifiers: strings requesting specific output
%       spbCenDistance : plots histogram, shows mean/std/median/%out of
%           distance of cen to closest spb.
%           output: movieIdx, t, cenIdx, distance, distanceSigma
%
% OUTPUT idlistList with additional fields named as input
%
% REMARKS
%
% created with MATLAB ver.: 7.6.0.324 (R2008a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 08-Apr-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get input

if nargin > 0 && isstruct(varargin{1})
    % it's an idlistList
    idlistList = varargin{1};
    varargin(1) = [];
else

    % load data
    idlistList = loadIdlistList(cdBiodata(4),...
        struct('check','ask'));
end

nIdlists = length(idlistList);
if nIdlists == 0
    return
end


%% General Info

% read distances, tagOrder etc.
for i=1:nIdlists

    [tagExists,tagOrder] = ismember({'spb1','cen1','spb2','cen2'},idlistList(i).idlist(1).stats.labelcolor);

    % calculate tagOrder for distance calculation
    idTagOrder = tagOrder(tagOrder>0);
    idTagOrder = idTagOrder - min(idTagOrder) + 1;
    nTags = length(idTagOrder);
    idlistList(i).nTags = nTags;
    
    idl = idlistList(i).idlist;
    idlistList(i).goodTimes = catStruct(1,'idl.linklist(1)');

    if checkIdlist(idlistList(i).idlist,1)
        % new idlist
        [idlistList(i).distance, idlistList(i).distanceUnitVectors] = ...
            idlist2distMat(idlistList(i).idlist, idlistList(i).dataProperties,[],[],idTagOrder);
    else

        % calculate distances, distanceVectors without idlist2distMat - we need to
        % be able to use old idlists. Since we don't worry about uncertainties,
        % it's still fairly straightforward
        linklists = cat(3,idlistList(i).idlist.linklist);
        linklists = linklists(idTagOrder,:,:);
        goodTimes = squeeze(linklists(1,1,:));
        nTimepoints = length(idlistList(i).idlist);

        idlistList(i).distance = repmat(NaN,[nTags,nTags,nTimepoints]);
        idlistList(i).distanceUnitVectors = repmat(NaN,[nTags,nTags,nTimepoints,3]);
        for t = goodTimes'
            % get distance matrix, distanceVectorMatrix. Divide
            % distanceVectorMatrix by distance to get normed vectors
            % use the same ordering as idlist2distMat
            [idlistList(i).distance(:,:,t),distanceVectorMatrix] =...
                distMat(linklists(:,9:11,linklists(1,1,:)==t));
            idlistList(i).distanceUnitVectors(:,:,t,:) = ...
                permute(...
                distanceVectorMatrix./repmat(idlistList(i).distance(:,:,t),[1,1,3]),[1,2,4,3]);
        end
    end
end


%% Specific Info

if nargin == 0
    return
end

if any(strmatch('spbCenDistance',varargin))
    for i=1:nIdlists
        % fill idlist.spbCenDistance
        %

        % switch nTags
        switch idlistList(i).nTags
            case 2
                % easy
                dc1 = reshape(idlistList(i).distance(2,1,:),[],1,1);
                sc1 = reshape(idlistList(i).distance(1,2,:),[],1,1);
                ntp = length(dc1);
                idlistList(i).spbCenDistance = [ones(ntp,1)*i,(1:ntp)',...
                    ones(ntp,1)*1,dc1,sc1];
            case 3
                % have to compare distances.
                dc11 = reshape(idlistList(i).distance(2,1,:),[],1,1);
                sc11 = reshape(idlistList(i).distance(1,2,:),[],1,1);
                dc12 = reshape(idlistList(i).distance(3,2,:),[],1,1);
                sc12 = reshape(idlistList(i).distance(2,3,:),[],1,1);

                dc1 = min(dc11,dc12);
                ntp = length(dc1);
                sc1=sc11;
                sc1(dc12==dc1) = sc12(dc12==dc1);

                idlistList(i).spbCenDistance = [ones(ntp,1)*i,(1:ntp)',...
                    ones(ntp,1)*1,dc1,sc1];
            case 4
                % do the same as case 3, but for two cen tags
                dc11 = reshape(idlistList(i).distance(2,1,:),[],1,1);
                sc11 = reshape(idlistList(i).distance(1,2,:),[],1,1);
                dc12 = reshape(idlistList(i).distance(3,2,:),[],1,1);
                sc12 = reshape(idlistList(i).distance(2,3,:),[],1,1);

                dc1 = min(dc11,dc12);
                ntp = length(dc1);
                sc1=sc11;
                sc1(dc12==dc1) = sc12(dc12==dc1);

                dc21 = reshape(idlistList(i).distance(4,1,:),[],1,1);
                sc21 = reshape(idlistList(i).distance(1,4,:),[],1,1);
                dc22 = reshape(idlistList(i).distance(4,3,:),[],1,1);
                sc22 = reshape(idlistList(i).distance(3,4,:),[],1,1);

                dc2 = min(dc21,dc22);
                %ntp = length(dc2);
                sc2=sc21;
                sc2(dc22==dc2) = sc22(dc22==dc2);

                idlistList(i).spbCenDistance = [ones(ntp,1)*i,(1:ntp)',...
                    ones(ntp,1)*1,dc1,sc1;...
                    ones(ntp,1)*i,(1:ntp)',...
                    ones(ntp,1)*2,dc2,sc2;];
        end

        % do distance stats: robustMean, robustStd, outlier%, n
        [rm,rs,ii,oo] = robustMean(idlistList(i).spbCenDistance(:,4));
        oo = length(oo);ii = length(ii);
        idlistList(i).spbCenDistanceStats = [rm,rs,oo/(ii+oo),ii+oo];


    end

    allDist = cat(1,idlistList.spbCenDistance);
    allStats = cat(1,idlistList.spbCenDistanceStats);

    % display overall stats and histogram
    figure('Name','SPB-CEN distances')
    histogram(allDist(:,4),1,0);

    [rm,rs] = robustMean(allDist(:,4));
    fprintf('Average Robust Distance : %1.2f%s%1.2f um\n',rm,char(177),rs)
    [rm,rs] = robustMean(allStats(:,2));
    fprintf('Average Robust Sigma    : %1.2f%s%1.2f um\n\n\n',rm,char(177),rs)

    % individual stats
    for i=1:nIdlists
        fprintf('Stats for %s : %1.2f%s%1.2f %%out:%2.1f n=%i\n',...
            idlistList(i).name,idlistList(i).spbCenDistanceStats(1),...
            char(177),idlistList(i).spbCenDistanceStats(2:end));
    end

end
