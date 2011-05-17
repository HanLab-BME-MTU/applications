function array = addBorder(array,border,nanMask)
%ADDBORDER adds a mirrored border around 1-3d arrays
%
% SYNOPSIS: array = addBorder(array,border)
%
% INPUT array: d-dimensional array, where d can be 1-3
%		border : 1-by-d array indicating by how many pixels the array
%                should be extended. Usually, that is going to be half of
%                your array size.
%                Pass a 2-by-d array if there is unequal padding on the two
%                sides of the array
%       nanMask: (opt) logical array of the same size as the image indicating
%                whether there are Nans. If there are any, the code 
%                adds a border around the masked image containing
%                the values (min+max)/2. Default: []
%
% OUTPUT array: array with added borders
%
% REMARKS addBorder reflects the array to pad it
%
% created with MATLAB ver.: 7.6.0.324 (R2008a) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 26-May-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if nargin < 2 || isempty(array) || isempty(border)
    error('please supply at least two non-empty input arguments')
end

nDims = ndims(array);
if nDims == 1
    array = array(:);
end
if nDims > 3
    error('please supply maximum 3D arrays')
end
if size(border,2) ~= nDims
    if size(border,2) > nDims && all(all(border(:,nDims+1:end))) == 0
        border = border(:,1:nDims);
    else
        error('border needs a column for every dimension')
    end
end
if size(border,1) == 1
    border = [border;border];
end
if nargin < 3 || isempty(nanMask)
    nanMask = [];
else
    if ~any(any(any(nanMask,1),2),3)
        nanMask = [];
    end
end

% do nanMask-stuff before flipping
if ~isempty(nanMask)
    % add border in nan-masked image. Do this per slice for speed (and because
    % there is no need for 3d right now. Throw a warning if there could be a
    % problem

    if nDims > 2 && any(all(all(nanMask,1),2))
        warning('ADDBORDER:NANSUPPORT','3D-masks are not fully supported yet')
    end
    % for slices that have NaNs, dilate/erode and take the average
    nanList = find(any(any(nanMask,1),2));
    minArray = squeeze(min(min(array,[],1),[],2));
    maxArray = squeeze(max(max(array,[],1),[],2));
    % create structuring element based on border+1 (just to be safe)
    se = strel('rectangle',sum(border(:,1:min(nDims,2)),1)+1);
    for z = nanList'
        % dilate. Replace NaN with small value
        tmp1 = array(:,:,z);
        if ~all(isnan(tmp1(:))) % no 3D support test here
        tmp1(nanMask(:,:,z)) = minArray(z)-1;
        tmp1 = imdilate(tmp1,se);
        tmp1(tmp1==minArray(z)-1) = NaN; % put mask back on
        % erode. Replace NaN with large value
        tmp2 = array(:,:,z);
        tmp2(nanMask(:,:,z)) = maxArray(z)+1;
        tmp2 = imerode(tmp2,se);
        tmp2(tmp2==maxArray(z)) = NaN; % put mask back on

        % take average
        tmp1 = (tmp1+tmp2)/2;
        % put original data back into tmp
        tmp2 = array(:,:,z);
        tmp1(~nanMask(:,:,z)) = tmp2(~nanMask(:,:,z));
        % replace slice
        array(:,:,z) = tmp1;
        end
    end % loop z
end% nanMask
        
% flip along first dimension
array = cat(1,flipdim(array(1:border(1,1),:,:),1),...
    array,flipdim(array(end-border(2,1)+1:end,:,:),1));
if nDims > 1
    % flip along second dimension
    array = cat(2,flipdim(array(:,1:border(1,2),:),2),...
    array,flipdim(array(:,end-border(2,2)+1:end,:),2));
end

if nDims > 2
    % flip along third dimension
     array = cat(3,flipdim(array(:,:,1:border(1,3)),3),...
    array,flipdim(array(:,:,end-border(2,3)+1:end),3));
end

