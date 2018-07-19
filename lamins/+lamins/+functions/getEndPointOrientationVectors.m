function [ vectorOut ] = getEndPointOrientationVectors( BW, g )
%getEndPointOrientationVectors g

        import lamins.functions.*;

        nAngles = 36;

        endpts = find(bwmorph(BW,'endpoints'));

        nfilter = ones(3);
        nfilter(2,2) = 0;
        F = imfilter(double(BW),nfilter).*BW;

        cc = bwconncomp(F == 2 | F == 1);
        cc = orderPixelIdxList(cc);
        BW = labelmatrix(cc) > 0;
        F = F.*BW;

    aresp = shiftdim(g(2).a,2);


    candidates = findCandidateIntersections(cc.PixelIdxList,g(2).a);
    detectedIntersections = find(F > 2);
    % [DY,DX] = ind2sub(size(I),detectedIntersections);
    ucandidates = unique(vertcat(candidates{:},detectedIntersections));
    
    fits = fitVonMisesFischer2dparallel(aresp(:,ucandidates)',0.1,4);
    fits = sortVonMisesFischer2d(fits);
    % 
    % ratio = cellfun(@(x) x(4),fits)./cellfun(@(x) x(2),fits);
    % ratio_T = thresholdOtsu(ratio);

%     fittedAngles = false(nAngles,length(ucandidates));
%     for i=1:length(fits)
%         fittedAngles(mod(round((fits{i}(1:2:end-2)+pi/2)/pi*36),36)+1,i) = true;
%     end

    map(ucandidates) = 1:length(ucandidates);

    candidateAngles = false(nAngles,length(ucandidates));
    for i=1:nAngles
        candidateAngles(i,map(candidates{i})) = true;
    end

%     matchedAngles = fittedAngles & candidateAngles;
% matchedCandidates = any(matchedAngles);

% [MY,MX] = ind2sub(size(I),ucandidates(matchedCandidates));


% plotCandidateIntersections(candidates);

% endpts = find(bwmorph(BW,'endpoints'));
endpt_fits = fitVonMisesFischer2dparallel(aresp(:,endpts)',0.1,4);
endpt_fits = sortVonMisesFischer2d(endpt_fits);
% [EY,EX] = ind2sub(size(I),endpts);


% BWcc = bwconncomp(BW,8);

%% match

ucandidates = ucandidates(cellfun('length',fits) > 4);
fits = fits(cellfun('length',fits) > 4);

intersection_stats = getEndpointStatistics(ucandidates,cc,fits,Inf);
endpt_stats = getEndpointStatistics(endpts,cc,endpt_fits,1);

intersection_stats = [intersection_stats ; endpt_stats];

D = pdist2(endpt_stats,intersection_stats,@costEndpointIntersection);
[matchedCost,matchedIdx] = min(D,[],2);

validMatch = matchedCost < Inf;

vectorOut = [intersection_stats(matchedIdx(validMatch),2) - endpt_stats(validMatch,2) ... 
           intersection_stats(matchedIdx(validMatch),1) - endpt_stats(validMatch,1)];

end

