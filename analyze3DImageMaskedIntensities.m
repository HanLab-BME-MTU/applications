function bp = analyze3DImageMaskedIntensities(images,mask,skelGraph,showPlots)

%DOCUMENTATION NEEDEDDDDDDDDDAADADAJJDKHAKLJDHAKLJDHALJKHDAJKL I HATE
%SCIENCE!!!!!!!!!!!!!!

[M N P nChan] = size(images);

showPlots = true;
distBinSz = 1;

distX = bwdist(~mask);

distBins = 0:distBinSz:max(distX(:));
nDistBins = numel(distBins);

bp.wholeMaskMean = nan(nChan,nDistBins-1);
bp.wholeMaskSTD = nan(nChan,nDistBins-1);

for k = 1:nChan
    
    for j = 1:numel(distBins)-1
        tmp = double(images(:,:,:,k));%Lazy, makes indexing easier
        bp.wholeMaskMean(k,j) = mean(tmp(distX > distBins(j) & distX < distBins(j+1) & mask));%Use mask in logical also for the first bin ... Also lazy...
        bp.wholeMaskSTD(k,j) = std(tmp(distX > distBins(j) & distX < distBins(j+1) & mask));
    end

end



if showPlots
    
    figure
    hold on
    
    chanCols = jet(nChan);
    chanStr = cell(nChan,1);
    for k = 1:nChan
        tmp = bp.wholeMaskMean(k,:) - nanmin(bp.wholeMaskMean(k,:));        
        tmp2 = bp.wholeMaskSTD(k,:) / nanmax(tmp);
        tmp = tmp ./ nanmax(tmp);
        plot(distBins(1:end-1),tmp,'color',chanCols(k,:));        
        chanStr{k} = ['Channel ' num2str(k)];
    end
    legend(chanStr{:})
    %SOOOO LAZY AND TIRED RIGHT NOW - DO IT THIS WAY SO LEGEND MAKES SENSE
    for k = 1:nChan
        tmp = bp.wholeMaskMean(k,:) - nanmin(bp.wholeMaskMean(k,:));        
        tmp2 = bp.wholeMaskSTD(k,:) / nanmax(tmp);
        tmp = tmp ./ nanmax(tmp);        
        plotTransparent(distBins(1:end-1),tmp,tmp2,chanCols(k,:),.3,0);        
    end
    
    
    xlabel('Distance from membrane, pixels')
    ylabel('Normalized Average Fluoresence')
    title('Distance - Intensity Profile, Whole-Cell')
    
    
end
    
