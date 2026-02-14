function saveRepresentativeSegmentationImages(md, outStruct, saveDir, prefix)

if ~exist(saveDir,'dir'); mkdir(saveDir); end

D = outStruct.images.DAPIproj;
R = outStruct.images.NFKBproj;

nucLabel  = outStruct.masks.nucLabel;
cytoLabel = outStruct.masks.cytoLabel;

% normalize to uint8
D8 = im2uint8(mat2gray(D));
R8 = im2uint8(mat2gray(R));

imwrite(D8, fullfile(saveDir, sprintf('%s_DAPI.png', prefix)));
imwrite(R8, fullfile(saveDir, sprintf('%s_NFkB.png', prefix)));

% outlines
nucB  = bwperim(nucLabel>0);
cytoB = bwperim(cytoLabel>0);

rgb = repmat(D8,1,1,3);

% nucleus outline: white
rgb(:,:,1) = max(rgb(:,:,1), uint8(255)*uint8(nucB));
rgb(:,:,2) = max(rgb(:,:,2), uint8(255)*uint8(nucB));
rgb(:,:,3) = max(rgb(:,:,3), uint8(255)*uint8(nucB));

% cyto outline: yellow
rgb(:,:,1) = max(rgb(:,:,1), uint8(255)*uint8(cytoB));
rgb(:,:,2) = max(rgb(:,:,2), uint8(255)*uint8(cytoB));

imwrite(rgb, fullfile(saveDir, sprintf('%s_DAPI_overlay.png', prefix)));

% outline only
ol = zeros(size(D8,1), size(D8,2), 3, 'uint8');
ol(:,:,1) = uint8(255)*uint8(cytoB);
ol(:,:,2) = uint8(255)*uint8(cytoB);
ol(:,:,1) = max(ol(:,:,1), uint8(255)*uint8(nucB));
ol(:,:,2) = max(ol(:,:,2), uint8(255)*uint8(nucB));
ol(:,:,3) = max(ol(:,:,3), uint8(255)*uint8(nucB));

imwrite(ol, fullfile(saveDir, sprintf('%s_outlineOnly.png', prefix)));

% meta
try
    fid = fopen(fullfile(saveDir, sprintf('%s_meta.txt', prefix)), 'w');
    fprintf(fid, 'MD path: %s\n', md.movieDataPath_);
    if isfield(outStruct,'params_px')
        fprintf(fid, 'umPerPx: %.6f\n', outStruct.params_px.umPerPx);
        fprintf(fid, 'CytoInnerPx: %d\n', outStruct.params_px.CytoInnerRadiusPx);
        fprintf(fid, 'CytoOuterPx: %d\n', outStruct.params_px.CytoOuterRadiusPx);
        fprintf(fid, 'BgOuterPx: %d\n', outStruct.params_px.BgOuterRadiusPx);
    end
    fclose(fid);
end
end