%% Load in data
load('project/140515/original reconstructed images/ag_080712wt_Reconstructed 3/ag_080712wt_Reconstructed 3.mat');
cellreader = CellReader(CachedReader(MD));
% [res, theta, nms] = steerableDetector(double(cellreader{3,1,15}),2,1);

%% Do steerableFilter detection on 3rd channel (Lamin B1)
nms = cell(1,20);
for i=1:20;
    [res, theta, nms{i}] = steerableDetector(double(cellreader{3,1,i}),2,1);
end;
I = cellreader;
% do contrast adjustment on the pictures
J = cellfun(@imadjust,I(1:4,1,15),'UniformOutput',0);
J = [J{:} {nms{15} > 0.5}];
K = im2double(cat(3,J{1:5}));
L = shiftdim(K(:,:,1:4),2);
L(:,nms{15} > 0.5) = 0;
K = shiftdim(L,1);
% imwrite(reshape(reshape(K(:,:,1:4),[],4)/1.3*jet(4),[1024 1024 3]),'frame15comp.png');

%% Find Branch points
bps = bwmorph(bwmorph(nms{15} > 0.5,'skel'),'branchpoints');
[r,c] = find(bps);

%% Histogram of pairwise distance
% dist{15} 


%% Plot
figure;
imshow(reshape(reshape(K(:,:,1:4),[],4)/1.3*jet(4),[1024 1024 3]));
hold on;
plot(c,r,'ro','MarkerFaceColor','r');

%% Scan parameters
sigma = 1:0.5:10.5;
scan = cell(1,length(sigma));
Q = cellreader{3,1,15};
for s=1:length(sigma);
    [~,~,scan{s}] = steerableDetector(double(Q),2,sigma(s));
end;

%% Show Parameter scan
colorStack(scan)

%% Show Z color stack
colorStack(I(3,1,:),'rgb',true,'imadjust',0)

%% Steerable filter on all lamin channels
goodsigma = sigma(12);
nms = cell(1,getSizeC(cellreader));
parfor c=1:getSizeC(cellreader)
    nms{c} = cell(1,getSizeT(cellreader));
    for t=1:getSizeT(cellreader)
        nms{c}{t} = cell(1,getSizeZ(cellreader));
        for z=1:getSizeZ(cellreader)
            [~,~,nms{c}{t}{z}] = steerableDetector(double(cellreader{c,t,z}),2,goodsigma);
        end
    end
end