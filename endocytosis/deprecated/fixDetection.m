% Francois Aguet, 08/05/2011

function fixDetection(data)

% fields to repair/copy
fnames = {'x', 'y', 'A', 'c',...
    'x_pstd', 'y_pstd', 'A_pstd', 'c_pstd',...
    'sigma_r', 'SE_sigma_r', 'RSS', 'pval_Ar'};

for k = 1:length(data)
    frameInfo = load([data(k).source 'Detection' filesep 'detection_v2.mat']);
    frameInfo = frameInfo.frameInfo;
    mCh = strcmp(data(k).channels, data(k).source);
    nCh = length(data(k).channels);
    sigma = getGaussianPSFsigma(data(k).NA, data(k).M, data(k).pixelSize, data(k).markers{mCh});
    for f = 1:data(k).movieLength
        for c = setdiff(1:nCh,mCh)
            nanIdx = find(isnan(frameInfo(f).x(c,:)));
            if ~isempty(nanIdx)
                img = double(imread(data(k).framePaths{c}{f}));
                
                xi = frameInfo(f).x_init(nanIdx);
                yi = frameInfo(f).y_init(nanIdx);
                
                pstruct = fitGaussians2D(img, xi, yi, [], sigma*ones(1,numel(nanIdx)), [], 'Ac');
                pstructLoc = fitGaussians2D(img, xi, yi, [], sigma*ones(1,numel(nanIdx)), [], 'xyAc');
                idx = sqrt((xi-pstructLoc.x).^2 + (yi-pstructLoc.y).^2) < 3*sigma & pstructLoc.A > pstruct.A;
                
                if sum(idx)~=0 % valid localization
                    for n = 1:length(fnames)
                        frameInfo(f).(fnames{n})(c,nanIdx(idx)) = pstructLoc.(fnames{n})(idx);
                    end
                end
                if sum(idx)~=length(nanIdx)
                    for n = 1:length(fnames)
                        frameInfo(f).(fnames{n})(c,nanIdx(~idx)) = pstruct.(fnames{n})(~idx);
                    end
                end
            end
        end
    end
    save([data(k).source 'Detection' filesep 'detection_v2.mat'], 'frameInfo');
end
