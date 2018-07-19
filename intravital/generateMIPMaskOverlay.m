function [imResult] = generateMIPMaskOverlay(imInput, masks, maskColors, maskAlphas, varargin)

    p = inputParser;
    p.CaseSensitive = false;
    p.addRequired('imInput', @(x) (isnumeric(x) && ndims(x) == 3) );
    p.addRequired('masks', @(x) ((all(size(x) == size(imInput))) || iscell(x)));
    p.addRequired('maskColors', @(x) (isnumeric(x)));
    p.addRequired('maskAlphas', @(x) (isnumeric(x)));
    p.addParamValue('spacing', ones( 1, ndims(imInput) ), @(x) (isnumeric(x) && ismember(numel(x), [1,ndims(imInput)])) );
    p.addParamValue('gapBetweenMIPs', 10, @(x) (isscalar(x) && x > 0));
    p.parse(imInput, masks, maskColors, maskAlphas, varargin{:});

    PARAMETERS = p.Results;
    
    if ~iscell( masks )
        masks = { masks };
    end    
    
    gapSize = PARAMETERS.gapBetweenMIPs;
    imStackMIP = cell(1,3);
    
    for dim = 1:3

        imCurMIP = mat2gray(squeeze(max(imInput, [], dim)));
        
        if dim < 3
            outRows = size(imCurMIP, 1);
            outCols = round( size(imCurMIP, 2) * PARAMETERS.spacing(3) / PARAMETERS.spacing(dim) );
            imCurMIP = imresize(imCurMIP, [outRows, outCols]);
        end
        
        if dim == 1
            imCurMIP = imCurMIP';
        end

        maskMIPs = cell( size(masks) );
        
        for i = 1:numel(masks)
            
            imCurMaskMIP = mat2gray(squeeze(max(masks{i}, [], dim)));
            
            if dim < 3
                outRows = size(imCurMaskMIP, 1);
                outCols = round( size(imCurMaskMIP, 2) * PARAMETERS.spacing(3) / PARAMETERS.spacing(dim) );
                imCurMaskMIP = imresize(imCurMaskMIP, [outRows, outCols]);
            end
            
            if dim == 1
                imCurMaskMIP = imCurMaskMIP';
            end
            
            maskMIPs{i} = imCurMaskMIP;
            
        end
        
        imStackMIP{dim} = genImageMaskOverlay(imCurMIP, maskMIPs, maskColors, maskAlphas);

    end
    
    % XY
    imResult = imStackMIP{3}; 

    % XZ
    imResult(end+1:end+gapSize, :, :) = 1;
    imResult = cat(1, imResult, imStackMIP{1});

    % YZ
    imResult(:,end+1:end+gapSize, :) = 1;
    imResult = cat(2, imResult, padarray(imStackMIP{2}, [size(imResult,1)-size(imStackMIP{2},1), 0, 0], 1.0, 'post'));

    % border all around
    imResult = padarray(imResult, [gapSize, gapSize, 0], 1.0);
        
end