function [imBlobMask] = generateBlobMask(imsize, blobLoc, blobRadii, spacing)

    seedPos = blobLoc * diag(spacing);
    kd = KDTreeSearcher(seedPos);

    pixelPos = ind2submat(size(imInput), (1:prod(imsize))') * diag(spacing);
    [closestSeedInd, distanceToSeed] = kd.knnsearch(pixelPos);

    imBlobMask = zeros(imsize);
    flagIsPixelInSeedVicinity = distanceToSeed <= blobRadii(closestSeedInd);
    imBlobMask( flagIsPixelInSeedVicinity ) = closestSeedInd( flagIsPixelInSeedVicinity );

end