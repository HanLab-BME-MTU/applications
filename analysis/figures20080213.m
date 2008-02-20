function varargout = figures20080213(figureId)
%FIGURES20080213 plots the figures discussed in the phone conference on Feb 13 2008
%
% SYNOPSIS: [data,idlistList] = figures20080213(figureId)
%
% INPUT figureId : id of figure to plot
%           1 - ratio of spb/cen intensities from G1 data
%           2 - relative position of cen for 3-spot movies
%
% OUTPUT 
%
% REMARKS
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 13-Feb-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if nargin == 0
    error('please input figureId')
end

% switch according to id
switch figureId
    case 1
        % ratio of spb/cen intensities from G1 data
        
        % load idlists
        % only take idlists/frames where spb1 and cen1 have been found and
        % aren't fused, and where there aren't more than two tags overall
        % currently, do not load loadStruct (and don't refit)
        done = false;
        idlistList = [];
        while ~done
            try
        tmp = loadIdlistList(cdBiodata(4),struct('check',{{7,'spb1','cen1';4,2,''}}),1,0);
        idlistList = [idlistList,tmp];
            catch
                done = true;
            end
        end
        
        
        
        % loop through list. Refit intensities, then read and collect
        nIdlists = length(idlistList);
        % s1c1int: movieIdx, t, distance, s1Int, c1Int, ratio
        % (intensities are stored in case I'd like to try exponential
        % fitting of intensities)
        s1c1int = zeros(200*nIdlists,6);
        scIdx = 0;
        
        for iIdlist = 1:nIdlists
            
            % refit intensities - don't do for now
            idlist = idlistList(iIdlist).idlist;
           
%             [idlist,dataProperties] = refitIntensities(idlistList(iIdlist).idlist,...
%                 idlistList(iIdlist).dataProperties,...
%                 fullfile(idlistList(iIdlist).loadStruct.moviePath,...
%                 idlistList(iIdlist).loadStruct.movieName),1,0);
            
            [tagExists,tagOrder] = ismember({'spb1','cen1'},idlist(1).stats.labelcolor);
            
            % read intensities
            linklists = cat(3,idlist.linklist);
            nTimepoints = size(linklists,3);
            dist = normList(squeeze(linklists(tagOrder(1),9:11,:)-linklists(tagOrder(2),9:11,:))');
            s1c1int(scIdx+1:scIdx+nTimepoints,2:5) = ...
                [squeeze(linklists(1,1,:)), dist,...
                squeeze(linklists(tagOrder(1),8,:)),...
                squeeze(linklists(tagOrder(2),8,:))];
            s1c1int(scIdx+1:scIdx+nTimepoints,1) = iIdlist;
            
            % update index
            scIdx = scIdx + nTimepoints;
        end
        
        % remove superfluous entries
        s1c1int(scIdx+1:end,:) = [];
        
        % take ratio
        s1c1int(:,6) = s1c1int(:,5) ./ s1c1int(:,4);
        
        % plot
        figure('Name','relative intensities')
        hold on
        for iIdlist = 1:nIdlists
            idx = find(s1c1int(:,1) == iIdlist);
            % get rid of multiple similar ratios from old tracked idlists
            try
                [dummy,idxIdx] = unique(s1c1int(idx,6),'first');
            catch
                 [dummy,idxIdx] = unique(s1c1int(idx,6));
            end
            plot(s1c1int(idx(idxIdx),3),s1c1int(idx(idxIdx),6),'.','Color',extendedColors(iIdlist));
        end
        xlabel('distance spb-cen')
        ylabel('ratio cenInt:spbInt')
        
        figure,histogram(unique(s1c1int(:,6)),1,1)
        
        varargout{1} = s1c1int;
           varargout{2} = idlistList;
    case 2
        
        % for metaphase 3-spot movies with spindles shorter than 1.4 um and no cen1*:
        % plot distribution of centromere distance
        done = false;
        idlistList = [];
        while ~done
            try 
                tmp = loadIdlistList(cdBiodata(4),struct('check',{{8,0,1.4;4,3,'';3,'',''}}),1,0);
 
        idlistList = [idlistList,tmp];
            catch
                done = true;
            end
        end
        
          
           
           % loop through idlistList. Read coords.
           nIdlists = length(idlistList);
           
           for i = 1:nIdlists
               % copied from cdMetaphaseGeometry
               [tagExists,tagOrder] = ismember({'spb1','cen1','cen2','spb2'},idlistList(i).idlist(1).stats.labelcolor);


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
               
               % distList: s1c1, s1s2
               idlistList(i).distList = [squeeze(idlistList(i).distance(2,1,:)),squeeze(idlistList(i).distance(3,1,:))];
               goodTimes = find(all(~isnan(idlistList(i).distList),2));
 % store index, time
               idlistList(i).i = i * ones(length(goodTimes),1);
               idlistList(i).t = goodTimes;
               
               % get scalar product of vectors for projection of s1c1 onto
               % s1s2
               s1c1vec = squeeze(idlistList(i).distanceUnitVectors(2,1,goodTimes,:));
               s1s2vec = squeeze(idlistList(i).distanceUnitVectors(3,1,goodTimes,:));
               proj = sum(s1c1vec.*s1s2vec,2);
               
               idlistList(i).distList = idlistList(i).distList(goodTimes,:);
               idlistList(i).distList(:,3) = (idlistList(i).distList(:,1).*proj)./idlistList(i).distList(:,2);
               if sum(idlistList(i).distList(:,3)>0.5) > 0.5 * length(goodTimes)
                   idlistList(i).distList(:,3) = 1-idlistList(i).distList(:,3);
               end
                   
           end
           
           % distance is interesting between 1 and 2
           s1c1dist = [cat(1,idlistList.i),cat(1,idlistList.t),cat(1,idlistList.distList)];
           
           figure('Name',sprintf('distance ratio for %i movies (%i tp)',nIdlists,size(s1c1dist,1))),
           histogram(s1c1dist(:,5),1,0)
ylabel('tp counts'),xlabel('s1c1:s1s2/s1s2 w/ flip')


           figure('Name',sprintf('projected distance ratio vs spindle length')),hold on
           for i=unique(s1c1dist(:,1)'),idx = s1c1dist(:,1) == i;plot(s1c1dist(idx,4),s1c1dist(idx,5),'.','Color',extendedColors(i)),end
           xlabel('s1s2 (\mum)'),ylabel('s1c1:s1s2/s1s2 w/ flip')
           figure('Name',sprintf('distance vs spindle length')),hold on
           for i=unique(s1c1dist(:,1)'),idx = s1c1dist(:,1) == i;plot(s1c1dist(idx,4),s1c1dist(idx,3),'.','Color',extendedColors(i)),end
           xlabel('s1s2 (\mum)'),ylabel('s1c1 (\mum)')
           figure('Name',sprintf('distance ratio vs spindle length')),hold on
           for i=unique(s1c1dist(:,1)'),idx = s1c1dist(:,1) == i;plot(s1c1dist(idx,4),s1c1dist(idx,3)./s1c1dist(idx,4),'.','Color',extendedColors(i)),end
           xlabel('s1s2 (\mum)'),ylabel('s1c1/s1s2')
           
           
           varargout{1} = s1c1dist;
           varargout{2} = idlistList;

    otherwise
        error('unrecognized figureId')
end