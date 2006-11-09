function idlistList = cdMetaphaseGeometry
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
            % cen2 is missing. Replace with cen1
            dTmp = idlistList(i).distance;
            dUvTmp = idlistList(i).distanceUnitVectors;
            % remove zero from tagOrder
%             tagOrder = tagOrder([1,2,4]); - no reordering!
tagOrder = 1:3;
            
            idlistList(i).distance = zeros(4,4,size(dTmp,3));
            idlistList(i).distance([1,2,4],[1,2,4],:) = dTmp(tagOrder,tagOrder,:);
            idlistList(i).distance(3,:,:) = idlistList(i).distance(2,:,:);
            idlistList(i).distance(:,3,:) = idlistList(i).distance(:,2,:);
            idlistList(i).distance(3,2,:) = 0;
            
            idlistList(i).distanceUnitVectors = zeros(size(dUvTmp) + [1,1,0,0]);
            idlistList(i).distanceUnitVectors([1,2,4],[1,2,4],:,:) = dUvTmp(tagOrder,tagOrder,:,:);
            idlistList(i).distanceUnitVectors(3,:,:,:) = idlistList(i).distanceUnitVectors(2,:,:,:);
            idlistList(i).distanceUnitVectors(:,3,:,:) = idlistList(i).distanceUnitVectors(:,2,:,:);
            idlistList(i).distanceUnitVectors(2,3,:,:) = 0;
            
            idlistList(i).nTags = 3;
        end
    else
        % remove data
        idlistList(i) = [];
    end
end

% count idlists again
nData = length(idlistList);

warning(s)

% loop through data, calculate geometry
for i=1:nData

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
    idlistList(i).meanOffsetVector = idlistList(i).meanOffsetVector;
    % relative projection of the centromere to centromere vector
    idlistList(i).centromereVectorProjection = dot(e_c1c2,e_s1s2,2);
    
    % spindle length
    idlistList(i).spindleLength = n_s1s2;
    
    % cenSpread: angle between cen-vector and spindle
    idlistList(i).maxCenSpreadAngle = ...
        max(acos(dot(e_s1c1,e_s1s2,2))*180/pi,acos(dot(e_s2c2,-e_s1s2,2))*180/pi);
    
end

% collect number of tags
nTags = cat(1,idlistList.nTags);

% full plotmatrix for 3, 4 and both

labels = {'dSPB (\mum)','dCen (\mum)','rel dCen','rel cencen-proj','mean off. (\mum)',...
    'rel mean off.', 'distanceRatio','max cenSpread angle', 'angle betw. off.-vec'};
labels3 = labels([1,5,6,7,8]);

% 3 spots
fh = figure('Name','PlotMatrix 3sp');
idx3 = nTags == 3;
dSpb = cat(1,idlistList(idx3).spindleLength);
mOff = cat(1,idlistList(idx3).meanOffsetVector);
mOffRel = cat(1,idlistList(idx3).n_meanOffsetVector);
dRat = cat(1,idlistList(idx3).distanceRatio);
aSpr = cat(1,idlistList(idx3).maxCenSpreadAngle);

[dummy,ah] = plotMatrix([dSpb,mOff,mOffRel,dRat,aSpr]);

for i=1:5
    set(get(ah(5,i),'xlabel'),'String',labels3{i});
    set(get(ah(i,1),'ylabel'),'String',labels3{i});
end


% 4 spots
fh = figure('Name','PlotMatrix4sp');
idx4 = nTags == 4;
dSpb = cat(1,idlistList(idx4).spindleLength);
dCen = cat(1,idlistList(idx4).centromereSeparation);
dCenRel = cat(1,idlistList(idx4).n_centromereSeparation);
dCenProj = cat(1,idlistList(idx4).centromereVectorProjection);
mOff = cat(1,idlistList(idx4).meanOffsetVector);
mOffRel = cat(1,idlistList(idx4).n_meanOffsetVector);
dRat = cat(1,idlistList(idx4).distanceRatio);
aSpr = cat(1,idlistList(idx4).maxCenSpreadAngle);
aOffV = cat(1,idlistList(idx4).offsetVectorAngle);


[dummy,ah] = plotMatrix([dSpb,dCen,dCenRel,dCenProj,mOff,mOffRel,dRat,aSpr,aOffV]);

for i=1:9
    set(get(ah(9,i),'xlabel'),'String',labels{i});
    set(get(ah(i,1),'ylabel'),'String',labels{i});
end











% plot
% fh = figure('Name','all');
% ylabels = {'angle between centromere offsets','normalized offset','path ratio'};
% for i=1:6
%     ah(i) = subplot(2,3,i);hold on
%     xlabel('normalized centromere separation')
%     ylabel(ylabels{mod(i-1,3)+1})
% end
% for i=1:nData
%     %figure('Name',idlistList(i).name)
%     %subplot(2,3,1)
%     axes(ah(1))
%     plot(idlistList(i).n_centromereSeparation,idlistList(i).offsetVectorAngle,'.','Color',extendedColors(i));
%     %subplot(2,3,2)
%     axes(ah(2))
%     plot(idlistList(i).n_centromereSeparation,idlistList(i).n_meanOffsetVector,'.','Color',extendedColors(i));
%     %subplot(2,3,3)
%     axes(ah(3))
%     plot(idlistList(i).n_centromereSeparation,idlistList(i).distanceRatio,'.','Color',extendedColors(i));
%     
%     %subplot(2,3,4)
%     axes(ah(4))
%     plot(diff(idlistList(i).n_centromereSeparation),diff(idlistList(i).offsetVectorAngle),'.','Color',extendedColors(i));
%     %subplot(2,3,5)
%     axes(ah(5))
%     plot(diff(idlistList(i).n_centromereSeparation),diff(idlistList(i).n_meanOffsetVector),'.','Color',extendedColors(i));
%     %subplot(2,3,6)
%     axes(ah(6))
%     plot(diff(idlistList(i).n_centromereSeparation),diff(idlistList(i).distanceRatio),'.','Color',extendedColors(i));
%     
%     
% end

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