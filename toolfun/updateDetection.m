function updateDetection(data)

nd = numel(data);

for k = 1:nd
    fprintf('Updating detection for %s ...', getShortPath(data(k)));
        
    % load data
    load([data(k).source 'Detection' filesep 'detection_v2.mat']);
    mCh = find(strcmp(data(k).source, data(k).channels));

    for f = 1:size(frameInfo)
        frameInfo(f).xCoord(:,2) = frameInfo(f).x_pstd(mCh,:);
        frameInfo(f).yCoord(:,2) = frameInfo(f).y_pstd(mCh,:);
        frameInfo(f).amp(:,2) = frameInfo(f).A_pstd(mCh,:);
    end
    save([data(k).source 'Detection' filesep 'detection_v2.mat'], 'frameInfo');
    fprintf(' done.\n');
end
