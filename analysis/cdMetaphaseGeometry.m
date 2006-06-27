function cdMetaphaseGeometry
%CDMETAPHASEGEOMETRY displays various metaphase geometry plots
%
% SYNOPSIS: cdMetaphaseGeometry
%
% INPUT
%
% OUTPUT
%
% REMARKS
%
% created with MATLAB ver.: 7.2.0.232 (R2006a) on Windows_NT
%
% created by: jdorn
% DATE: 22-Jun-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% read idlists, dataProperties. Only if 4 tags, no '?'.
idlistList = loadIdlistList(cdBiodata(4),...
    'length(idlist(1).stats.labelcolor) > 3 && isempty(strmatch(''?'',idlist(1).stats.labelcolor)) ');

nData = length(idlistList);

s  =  warning('query', 'all');
warning off MATLAB:divideByZero

% loop and get distances, vectors
for i=nData:-1:1
    % find order of tags
    [tagExists,tagOrder] = ismember({'spb1','cen1','cen2','spb2'},idlistList(i).idlist(1).stats.labelcolor);
    if ~all(tagExists)
        idlistList(i) = [];
    else
        if checkIdlist(idlistList(i).idlist,1)
            % new idlist
            [idlistList(i).distance, idlistList(i).distanceUnitVectors] = ...
                idlist2distMat(idlistList(i).idlist, idlistList(i).dataProperties);
        else

            % calculate distances, distanceVectors without idlist2distMat - we need to
            % be able to use old idlists. Since we don't worry about uncertainties,
            % it's still fairly straightforward
            linklists = cat(3,idlistList(i).idlist.linklist);
            goodTimes = squeeze(linklists(1,1,:));
            nTimepoints = length(idlistList(i).idlist);
            nTags = length(tagOrder);
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
        % re-sort distance, vectors
        idlistList(i).distance = idlistList(i).distance(tagOrder,tagOrder,:);
        idlistList(i).distanceUnitVectors = idlistList(i).distanceUnitVectors(tagOrder,tagOrder,:,:);
    end
end

% count idlists again
nData = length(idlistList);

warning(s)

% loop through data, calculate geometry
for i=1:nData
    % plot 2-by-3. Cols: angle between cen-tag-offset vectors, average
    % cen-tag-offset, ratio of s1-c1-c2-s1 to s1-s2.
    % Top row: 'static' data (cen-cen: elastic)
    % Bottom row: first difference (cen-cen: viscoelastic)

    % find cen-tag offset.
    
    % read vectors and distances
    e_s1s2 = squeeze(idlistList(i).distanceUnitVectors(1,4,:,:));
    e_s1c1 = squeeze(idlistList(i).distanceUnitVectors(1,2,:,:));
    e_s2c2 = squeeze(idlistList(i).distanceUnitVectors(4,3,:,:));
    e_c1c2 = squeeze(idlistList(i).distanceUnitVectors(2,3,:,:));
    n_s1c1 = squeeze(idlistList(i).distance(2,1,:));
    n_s2c2 = squeeze(idlistList(i).distance(4,3,:));
    n_c1c2 = squeeze(idlistList(i).distance(3,2,:));
    n_s1s2 = squeeze(idlistList(i).distance(4,1,:));
    v_s1c1 = e_s1c1 .* repmat(n_s1c1,1,3);
    v_s2c2 = e_s2c2 .* repmat(n_s2c2,1,3);

    % get cen1offsetVector
    v_cen1Offset = v_s1c1 - (e_s1s2 .* repmat(n_s1c1.*dot(e_s1s2,e_s1c1,2),1,3));
    [n_cen1Offset,e_cen1Offset] = normList(v_cen1Offset);
    
    % get cen2offsetVector
    v_cen2Offset = v_s2c2 - (e_s1s2 .* repmat(n_s2c2.*dot(e_s1s2,e_s2c2,2),1,3));
    [n_cen2Offset,e_cen2Offset] = normList(v_cen2Offset);
    
    % offsetVectorAngle is the angle between the two tag offset vectors
    idlistList(i).offsetVectorAngle = acos(dot(e_cen1Offset,e_cen2Offset,2))*180/pi;

    % meanOffsetVector is the average of the two offset lengths
    idlistList(i).meanOffsetVector = mean([n_cen1Offset,n_cen2Offset],2);
    
    % ratio of distance s1-c1-c2-s2 to s1-s2
    idlistList(i).distanceRatio = sum([n_s1c1,n_c1c2,n_s2c2],2)./n_s1s2;
    
    % centromere to centromere distance
    idlistList(i).centromereSeparation = n_c1c2;
    
    idlistList(i).n_centromereSeparation = n_c1c2./n_s1s2;
    idlistList(i).n_meanOffsetVector = idlistList(i).meanOffsetVector./n_s1s2;
    
    % relative projection of the centromere to centromere vector
    idlistList(i).centromereVectorProjection = dot(e_c1c2,e_s1s2,2);
end

% plot
fh = figure('Name','all');
ylabels = {'angle between centromere offsets','normalized offset','path ratio'};
for i=1:6
    ah(i) = subplot(2,3,i);hold on
    xlabel('normalized centromere separation')
    ylabel(ylabels{mod(i-1,3)+1})
end
for i=1:nData
    %figure('Name',idlistList(i).name)
    %subplot(2,3,1)
    axes(ah(1))
    plot(idlistList(i).n_centromereSeparation,idlistList(i).offsetVectorAngle,'.','Color',extendedColors(i));
    %subplot(2,3,2)
    axes(ah(2))
    plot(idlistList(i).n_centromereSeparation,idlistList(i).n_meanOffsetVector,'.','Color',extendedColors(i));
    %subplot(2,3,3)
    axes(ah(3))
    plot(idlistList(i).n_centromereSeparation,idlistList(i).distanceRatio,'.','Color',extendedColors(i));
    
    %subplot(2,3,4)
    axes(ah(4))
    plot(diff(idlistList(i).n_centromereSeparation),diff(idlistList(i).offsetVectorAngle),'.','Color',extendedColors(i));
    %subplot(2,3,5)
    axes(ah(5))
    plot(diff(idlistList(i).n_centromereSeparation),diff(idlistList(i).n_meanOffsetVector),'.','Color',extendedColors(i));
    %subplot(2,3,6)
    axes(ah(6))
    plot(diff(idlistList(i).n_centromereSeparation),diff(idlistList(i).distanceRatio),'.','Color',extendedColors(i));
    
    
end

% figure
% subplot(2,3,1)
%     plot(catStruct(1,'idlistList.n_centromereSeparation',NaN),catStruct(1,'idlistList.offsetVectorAngle',NaN),'.');
%     subplot(2,3,2)
%     plot(catStruct(1,'idlistList.n_centromereSeparation',NaN),catStruct(1,'idlistList.n_meanOffsetVector',NaN),'.');
%     subplot(2,3,3)
%     plot(catStruct(1,'idlistList.n_centromereSeparation',NaN),catStruct(1,'idlistList.distanceRatio',NaN),'.');
%     
%     subplot(2,3,4)
%     plot(diff(catStruct(1,'idlistList.n_centromereSeparation',NaN)),diff(catStruct(1,'idlistList.offsetVectorAngle',NaN)),'.');
%     subplot(2,3,5)
%     plot(diff(catStruct(1,'idlistList.n_centromereSeparation',NaN)),diff(catStruct(1,'idlistList.n_meanOffsetVector',NaN)),'.');
%     subplot(2,3,6)
%     plot(diff(catStruct(1,'idlistList.n_centromereSeparation',NaN)),diff(catStruct(1,'idlistList.distanceRatio',NaN)),'.');
%     