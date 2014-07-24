function [kmtLength,spindleLength,avgKmtLength,avgSpindleLength,idlistList] = cdPlotKinMtLength
%CDPLOTKINMTLENGTH plots kinetochore microtubule lengths
%
% SYNOPSIS: [kmtLength,spindleLength,avgKmtLength,avgSpindleLength,idlistList] = cdPlotKinMtLength
%
% INPUT
%
% OUTPUT kmtLength  : kinetochore microtubule length in um
%        spindleLength : corresponding spindle length in um
%        avg...     : averages of kmt (1,2) and spindle length per movie
%        idlistList : structure containing lots of information. The movie
%                     names corresponding to the numbers in the scatter
%                     plot are idlistList(#).name
%
% REMARKS kmtLength and spindleLength will contain NaN for wherever there
%         was no data
%         The function will warn automatically of movies with average kMT
%         length > 0.8 um
%
% created with MATLAB ver.: 7.7.0.2162 (R2008b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 22-Jul-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read idlists, dataProperties. Only if 4 tags, no '?'.
idlistList = loadIdlistList(cdBiodata(4),...
    'length(idlist(1).stats.labelcolor) > 2 && isempty(strmatch(''?'',idlist(1).stats.labelcolor)) ');

nData = length(idlistList);

s  =  warning('query', 'all');
warning off MATLAB:divideByZero

% loop and get distances, vectors
for i=nData:-1:1
    % find order of tags
    [tagExists,tagOrder] = ismember({'spb1','cen1','cen2','spb2'},idlistList(i).idlist(1).stats.labelcolor);
    if all(tagExists([1,2,4]))
        
        % calculate tagOrder for distance calculation
        idTagOrder = tagOrder(tagOrder>0);
        idTagOrder = idTagOrder - min(idTagOrder) + 1;
        
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
            nTags = length(idTagOrder);
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
        % re-sort distance, vectors. Check whether there may be just three
        % tags
        if tagExists(3)
            % all tags, no problem - removed reordering
            %         idlistList(i).distance = idlistList(i).distance(tagOrder,tagOrder,:);
            %         idlistList(i).distanceUnitVectors = idlistList(i).distanceUnitVectors(tagOrder,tagOrder,:,:);
            
            idlistList(i).nTags = 4;
        else
            % cen2 is missing. Replace with NaN
            
            dTmp = idlistList(i).distance;
            dUvTmp = idlistList(i).distanceUnitVectors;
            % remove zero from tagOrder
            %             tagOrder = tagOrder([1,2,4]); - no reordering!
            tagOrder = 1:3;
            
            idlistList(i).distance = zeros(4,4,size(dTmp,3));
            idlistList(i).distance([1,2,4],[1,2,4],:) = dTmp(tagOrder,tagOrder,:);
            idlistList(i).distance(3,:,:) = NaN;
            idlistList(i).distance(:,3,:) = NaN;
            
            idlistList(i).distanceUnitVectors = zeros(size(dUvTmp) + [1,1,0,0]);
            idlistList(i).distanceUnitVectors([1,2,4],[1,2,4],:,:) = dUvTmp(tagOrder,tagOrder,:,:);
            idlistList(i).distanceUnitVectors(3,:,:,:) = NaN;
            idlistList(i).distanceUnitVectors(:,3,:,:) = NaN;
            
            idlistList(i).nTags = 3;
        end
    else
        % remove data
        idlistList(i) = [];
    end
end

% count idlists again
nData = length(idlistList);

% collect distances
spindleLength = repmat(squeeze(catStruct(3,'idlistList.distance(4,1,:)')),2,1);
kmtLength = squeeze(cat(3,catStruct(3,'idlistList.distance(2,1,:)'),...
    catStruct(3,'idlistList.distance(4,3,:)')));

% plot
figure,scattercloud(spindleLength,kmtLength,[],[],'.k')
ylim([0,1.05*max(kmtLength)])
% write movie number
avgSpindleLength = zeros(nData,1);
avgKmtLength = zeros(nData,2);
for i=1:nData,
    avgSpindleLength(i)=nanmean(idlistList(i).distance(4,1,:));
    avgKmtLength(i,1)=nanmean(idlistList(i).distance(2,1,:));
    avgKmtLength(i,2)=nanmean(idlistList(i).distance(4,3,:));
end
hold on
for i=1:nData
    text(avgSpindleLength(i),avgKmtLength(i,1),num2str(i+0.1),'Color','r','FontWeight','bold')
    text(avgSpindleLength(i),avgKmtLength(i,2),num2str(i+0.2),'Color','r','FontWeight','bold')
end
% list movies with avg>0.8um
[u,v] = find(avgKmtLength>0.8);
if ~isempty(u)
    fprintf('movies to check (avg spb-cen distance > 0.8 um\n')
    for i=1:length(u)
        fprintf('%s cen%i \n',idlistList(u(i)).name,v(i))
    end
end


% clean up
warning(s)
