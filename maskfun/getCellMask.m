%[mask] = getCellMask(data, varargin)
%
% Notes: independent of the masking strategy, the detection must be run before

% Francois Aguet, 02/17/2012

function [mask] = getCellMask(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('Channel', 1, @isposint);
ip.addParamValue('Connect', true, @islogical);
ip.addParamValue('Display', false, @islogical);
ip.addParamValue('ShowHistogram', false, @islogical);
ip.addParamValue('ModeRatio', 0.8, @isscalar);
ip.addParamValue('Mode', 'intensity', @(x) any(strcmpi(x, {'intensity', 'density'})));
ip.parse(data, varargin{:});

nd = numel(data);
mask = cell(1,nd);
ch = ip.Results.Channel;

se = strel('disk', 5);

for i = 1:nd
    maskPath = [data(i).source 'Detection' filesep 'cellmask.tif'];
    aipPath = [data(i).source 'Detection' filesep 'avgProj.mat'];
    
    if ~(exist(maskPath, 'file')==2) || ip.Results.Overwrite
        
        switch ip.Results.Mode
            case 'intensity'
                % load max. 100 frames
                frameRange = unique(round(linspace(1, data(i).movieLength, 100)));
                aip = zeros(data(i).imagesize);
                mproj = zeros(data(i).imagesize);
                
                if iscell(data(i).framePaths{ch})
                    parfor f = 1:numel(frameRange)
                        frame = double(imread(data(i).framePaths{ch}{frameRange(f)})); %#ok<PFBNS>
                        dmask = 0~=double(imread(data(i).maskPaths{frameRange(f)}));
                        
                        dmask = imdilate(dmask, se);
                        frame(dmask) = 0;
                        aip = aip + frame;
                        mproj = mproj + dmask;
                    end
                else
                    info = imfinfo(data(i).framePaths{ch});
                    minfo = imfinfo(data(i).maskPaths);
                    parfor f = 1:numel(frameRange)
                        frame = double(readtiff(data(i).framePaths{ch}, frameRange(f), info)); %#ok<PFBNS>
                        dmask = 0~=double(readtiff(data(i).maskPaths, frameRange(f), minfo));
                        dmask = imdilate(dmask, se);
                        frame(dmask) = 0;
                        aip = aip + frame;
                        mproj = mproj + dmask;
                    end
                end
                
                aip = aip./(numel(frameRange)-mproj);
                % fill to enable filtering
                aip(isnan(aip)) = prctile(aip(:), 95);
                save(aipPath, 'aip');
                
                mask{i} = maskFromFirstMode(aip, 'Connect', ip.Results.Connect,...
                    'Display', ip.Results.ShowHistogram, 'ModeRatio', ip.Results.ModeRatio);
                
                mask{i} = imfill(mask{i}, 'holes');
                
                imwrite(uint8(mask{i}), maskPath, 'tif', 'compression' , 'lzw');
            case 'density'
                
    
        end    
    else
        mask{i} = double(imread(maskPath));
    end
end

if ip.Results.Display
    for i = 1:nd
        if ~isempty(mask{i})
            [ny,nx] = size(mask{i});
            B = bwboundaries(mask{i});
            B = cellfun(@(i) sub2ind([ny nx], i(:,1), i(:,2)), B, 'unif', 0);
            B = cell2mat(B);
            bmask = zeros([ny nx]);
            bmask(B) = 1;
            bmask = bwmorph(bmask, 'dilate');
            
            aipPath = [data(i).source 'Detection' filesep 'avgProj.mat'];
            load(aipPath);
            aip = scaleContrast(aip);
            aip(bmask==1) = 0;
            overlay = aip;
            overlay(bmask==1) = 255;
            overlay = uint8(cat(3, overlay, aip, aip));
            figure; imagesc(overlay); axis image; colormap(gray(256)); colorbar;
        end
    end
end

if nd==1
    mask = mask{1};
end
