function minD = distance2edge(path_to_spk1, path_to_spk2, path_to_mask, pixel_size)
% input:
%
% path_to_spk1: path to a fsmCenter project directory where the first
% speckle set is located.
%
% path_to_spk2: path to a fsmCenter project directory where the second
% speckle set is located.
%
% path_to_mask: path to a fsmCenter project directory where masks should
% be used.
%
% pixel_size: (in nanometers)
%
% output:
%
% minD: a matrix of size nFrames x 2 where nFrames is the number
% of frames. minD(i, 1) (resp. minD(i, 2)) corresponds to the min
% distance from speckles of set 1 (resp. set 2) to the edges(i).

N = 50;

files_spk1 = dir([path_to_spk1, filesep, 'cands*']);
files_spk2 = dir([path_to_spk2, filesep, 'cands*']);
files_mask = dir([path_to_mask, filesep, '*.tif']);

if (size(files_spk1, 1) ~= size(files_spk2, 1) || size(files_spk1, 1) ~= size(files_mask, 1))
    error('number of frames in subdirectories are not the same.');
end

nFrames = size(files_spk1, 1);

minD = zeros(nFrames, 2);

for frame=1:nFrames
    disp(frame);
    
    % Read the mask image
    BW = imread([path_to_mask, filesep, files_mask(frame).name]);
    
    % Read the speckle sets
    S1 = load([path_to_spk1, filesep, files_spk1(frame).name]);
    S2 = load([path_to_spk2, filesep, files_spk2(frame).name]);
    
    p1 = zeros(numel(S1.cands), 2);
    p2 = zeros(numel(S2.cands), 2);

    for i = 1:numel(S1.cands)
        p1(i, :) = S1.cands(i).Lmax;
    end
    
    for i = 1:numel(S2.cands)
        p2(i, :) = S2.cands(i).Lmax;
    end
    
    % Compute the distance transform   
    D = bwdist(max(BW(:)) - BW) * pixel_size;

    % Compute linear indices of speckles
    idxS1 = sub2ind(size(D), p1(:, 1), p1(:, 2));
    idxS2 = sub2ind(size(D), p2(:, 1), p2(:, 2));
    
    % Sort speckle according to their distance to the front
    minD1 = sort(D(idxS1));
    minD2 = sort(D(idxS2));
    
    % Compute the mean distance to the front of the first N speckles
    % closest to the front.

    idxNNZ1 = find(minD1, 1, 'first');
    idxNNZ2 = find(minD2, 1, 'first');
    
    minD(frame, 1) = mean(minD1(idxNNZ1:N));
    minD(frame, 2) = mean(minD2(idxNNZ2:N));
end

figure, plot(minD);
set(gca,'XTick',1:nFrames);
xlabel('Frames');
ylabel('Distance (nm)');
legend('Set 1', 'Set 2');
title('Minimum speckles distance to the edge');
end