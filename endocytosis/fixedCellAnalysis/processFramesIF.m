%processFramesIF(data, varargin) runs CCP detection on all channels of fixed-cell
% data sets for colocalization studies. To determine the level of chance
% colocalization, a set of Gaussian fits at random locations is also calculated.
%
% Input:
%    data : structure with fields   
%          .channels : cell array of paths to TIFF files
%          .results  : path to results file
%
% Options:
%    'Overwrite' : true|{false}
%
% The output structure written to the results file contains the fields
%  .ps     : result of pointSourceDetection at CCP locations, for all channels
%  .psRand : result of Gaussian fitting at random locations, for all channels
%  .mask   : cell mask

% Francois Aguet, 01/2014

function processFramesIF(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addParameter('Overwrite', false, @islogical);
ip.parse(data, varargin{:});

nd = numel(data);

for i = 1:nd
    
    if ~(exist(data(i).results, 'file')==2) || ip.Results.Overwrite
        fprintf('Processing %s ... ', getDirFromPath(data(i).results));
        
        nc = numel(data(i).channels);
        
        sigma = 1.5+zeros(1,nc); % change to input
        
        ch = cell(1,nc);
        for c = 1:nc
            ch{c} = imread(data(i).channels{c});
        end
        
        mask = maskFromFirstMode(ch{1});
        CC = bwconncomp(~mask, 8);
        np = cellfun(@numel, CC.PixelIdxList);
        mask = imfill(mask);
        bgIdx = vertcat(CC.PixelIdxList{np<100});
        mask(bgIdx) = 1;
        
        % run spot detection
        for c = 1:nc
            ps(c) = pointSourceDetection(ch{c}, sigma(c), 'Mode', 'xyAc', 'Mask', mask);
        end
        
        % random positions
        np = max(arrayfun(@(i) numel(i.x), ps));
        N = 10*np;
        w = max(ceil(4*sigma));
        % random positions
        xr = [];
        yr = [];
        [ny,nx] = size(mask);
        while numel(xr)<N
            xcand = 1+w+(nx-2*w-1)*rand(1,N);
            ycand = 1+w+(ny-2*w-1)*rand(1,N);
            idx = mask(sub2ind([ny nx], round(ycand), round(xcand)));
            xr = [xr xcand(idx~=0)]; %#ok<*AGROW>
            yr = [yr ycand(idx~=0)]; %#ok<*AGROW>
        end
        xr = xr(1:N);
        yr = yr(1:N);
        
        for c = 1:nc
            psRand(c) = fitGaussians2D(double(ch{c}), xr, yr, [], sigma(c), [], 'Ac'); %#ok<NASGU>
        end

        % save mask + detections
        save(data(i).results, 'ps', 'psRand', 'mask');
        fprintf('done.\n');
    end
end
