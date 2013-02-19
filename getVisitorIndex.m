% TO DO: change to cell-specific calculation

function vidx = getVisitorIndex(lftData, mCh)

if nargin<2
    mCh = 1;
end

minLength = min(arrayfun(@(i) size(i.gapMat_Ia,2), lftData));

A = arrayfun(@(i) i.A(:,1:minLength,mCh), lftData, 'UniformOutput', false);
A = vertcat(A{:});
lftV = arrayfun(@(i) i.lifetime_s(i.catIdx==1), lftData, 'unif', 0);
lifetime_s_all = vertcat(lftV{:});

% startBufferA = arrayfun(@(i) i.sbA(:,:,mCh), lftData, 'UniformOutput', false);
% startBufferA = vertcat(startBufferA{:});

tx = 30;

% 95th percentile of the reference (above threshold) distribution
% pRef = prctile(A(lifetime_s_all>=tx,1:3), 95, 1);

% pAbove = prctile(A(lifetime_s_all>=tx,1:3), 95, 1);
% pBelow = prctile(intMat_Ia_all(lifetime_s_all>=tx,1:3), 2.5, 1);

% pBuffer = prctile(startBufferA(lifetime_s_all>=tx,:), 5);


% 95th percentile of the reference (above threshold) distribution
pRef = prctile(A(lifetime_s_all>=tx,1:3,mCh), 95, 1);

nd = numel(lftData);
vidx = cell(1,nd);
for i = 1:nd
    vidx{i} = sum(lftData(i).A(:,1:3,mCh)>repmat(pRef, [size(lftData(i).A,1) 1]),2)>0 & lftV{i}<tx;
    
    % median of first three frames larger than median of last three frames
    %medStart = median(lftData(i).A(:,1:5),2);
    %medEnd = arrayfun(@(k) median(lftData(i).A(k,lftData(i).trackLengths(k)+(-2:0))), 1:size(lftData(i).A,1))';
    
    %maxStart = nanmax(lftData(i).A(:,1:5,mCh),[],2);
    %maxEnd = arrayfun(@(k) max(lftData(i).A(k,lftData(i).trackLengths(k)+(-4:0),mCh)), 1:size(lftData(i).A,1))';
    %maxEnd = nanmax(lftData(i).A(:,6:end),[],2);
    
    %diffIdx = maxStart>maxEnd;
    %idxLft = idxLft & diffIdx';
    
    %idxLft = (lftData(i).sbA(:,3)<pBuffer(3) | lftData(i).sbA(:,4)<pBuffer(4) | lftData(i).sbA(:,5)<pBuffer(5))' &...
    %    (lftData(i).A(:,1)>pAbove(1) | lftData(i).A(:,2)>pAbove(2) | lftData(i).A(:,3)>pAbove(3))';
    
    %idxLft = (lftData(i).sbA(:,4)<pBuffer(4) & lftData(i).sbA(:,5)<pBuffer(5))';
    
    %idxLft = sum(lftData(i).A(:,1:3)>repmat(pAbove, [size(lftData(i).A,1) 1]) |...
    %    lftData(i).A(:,1:3)<repmat(pBelow, [size(lftData(i).A,1) 1]),2)'>=1 & res(i).lft_all<tx;
    
    
    
end
