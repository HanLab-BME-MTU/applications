function W = computePairWeights()

% 6) For every pair of tracks, compute the probability of association

isValidPT = false(size(pairIdx,1),1);
weights = zeros(size(pairIdx,1),1);

% Compute x0, y0, amp0, sX0, sY0 and theta0, initial parameters set for each
% pair of track
% TODO

for iFrame = 1:nFrame
    % Read input image
    % TODO
    
    % Find which pair has any of its tracks living in iFrame
    indPT = find(pairIdx(:,1) == iFrame | pairIdx(:,2) == iFrame);
    
    for iiPT = 1:numel(indPT)
        iPT = indPT(iiPT);
        
        % Crop Image
        % TODO
        
        % Fit Model
        [params,stdParams,res] = fitAnisoGaussian2D(crop,[x0,y0,amp0,sX0,sY0,theta,C],'xyArstC');
        
        % Test Model
        % TODO
        % isValidPT(iPT) = true;
        
        % Assign BIC
        % TODO
        % weights(iPT) = BIC;
    end
end

% 7) Trim invalid pair of tracks

% Make additional tests here
% TODO
pairIdx = pairIdx(isValidPT,:);
weights = weights(isValidPT,:);

fprintf('Valid track pairs = %f %%\n',...
    nnz(isValidPT) * 100 / numel(hasOverlap));

% 8) Turn weights (BIC) into integer
% TODO
