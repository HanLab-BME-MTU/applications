function err = testTrackStackFlow(xShift,yShift,corLength,method,beadDensity,onlyImgs)
% scriptCompareCCWSvsSPC
% testTrackStackFlow finds error distribution of tracking algorithm in CDWS, CCWS
% and SPC for one-time-pass by making synthetic bead images
% example: err = testTrackStackFlow(0,0,17,'CDWS',5)
if nargin<6
    onlyImgs = false;
end
%% Bead characterization
NA = 1.4;
lambda = 647e-9;
camPixSize = 6.48e-6; % in m
objMag = 90;
sigma = getGaussianPSFsigma(NA,objMag,camPixSize,lambda); % this leads to 1.67
% sigma = 2.0;
%% Bead density characterization
xmax=100;
ymax=100;
pixSize = 72; % nm/pix 90x
% beadDensity = 5; % number per 10 pix x 10 pix
nPoints = xmax*ymax*beadDensity/100;
% bead_r = 50; % nm
%% Bead image creation
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
        'npoints', nPoints, 'Border', 'truncated','A',0.3+rand(1,nPoints),'pixelSize',pixSize);
refimg = refimg+0.05*rand(ymax,xmax)*max(refimg(:));
% figure, imshow(refimg,[])
%% Apply (sub-pixel) displacement - uniform
% nBeads = length(bead_x);
% xShift=0.1;
% yShift=0.1;
bead_ux = ones(size(bead_x))*xShift;
bead_uy = ones(size(bead_y))*yShift;

beadimg = simGaussianBeads(xmax,ymax, sigma,'x',bead_x+bead_ux,'y',bead_y+bead_uy,'A',Av, 'Border', 'truncated');
%% Noise addition (5%)
beadimg = beadimg+0.05*rand(ymax,xmax)*max(beadimg(:));
% figure, imshow(beadimg,[])
if onlyImgs
    imgPath = ['.' filesep 'imgPairU' num2str(xShift)];
    mkdir(imgPath)
    imgPath1 = [imgPath filesep 'img1.tif'];
    imgPath2 = [imgPath filesep 'img2.tif'];
    imwrite(refimg,imgPath1,'tif');
    imwrite(beadimg,imgPath2,'tif');
    return
end
    
%% beads sampling for tracking
% Select only beads which are min correlation length away from the border of the
% reference frame 
% corLength = 17;
beadsMask = true(size(refimg));
erosionDist = corLength;
beadsMask(erosionDist:end-erosionDist,erosionDist:end-erosionDist)=false;
indx=beadsMask(sub2ind(size(beadsMask),round(bead_y),round(bead_x)));
beads = [bead_x(~indx) bead_y(~indx)];
uBeads = [bead_ux(~indx) bead_uy(~indx)];
% hold on, plot(beads(:,1),beads(:,2),'ro')
%% Track with CDWS, CCWS and SPC
% track with CDWS
maxFlowSpeed = 5;
if strcmp(method,'FCTR')
    v = trackStackFlow(cat(3,refimg,beadimg),beads,...
        corLength,corLength,'maxSpd',maxFlowSpeed,...
        'mode','fast');
    err = ((v(:,1)-uBeads(:,1)).^2+(v(:,2)-uBeads(:,2)).^2).^0.5;
elseif strcmp(method,'CDWS')
    v = trackStackFlow(cat(3,refimg,beadimg),beads,...
        corLength,corLength,'maxSpd',maxFlowSpeed,...
        'mode','CDWS');
    err = ((v(:,1)-uBeads(:,1)).^2+(v(:,2)-uBeads(:,2)).^2).^0.5;
%% CCWS
elseif strcmp(method,'CCWS')
%     v = trackStackFlow(cat(3,refimg,beadimg),beads,...
%         corLength,corLength,'maxSpd',maxFlowSpeed,...
%         'mode','CCWS');
    pivPar = [];      % variable for settings
    pivData = [];     % variable for storing results
    im1 = refimg;
    im2 = beadimg;
    [pivPar, pivData] = pivParams(pivData,pivPar,'defaults');     
    [pivData] = pivAnalyzeImagePair(im1,im2,pivData,pivPar);
    v = [pivData.U(:) pivData.V(:)];
    nLen = length(v);
    err = ((v(:,1)-ones(nLen,1)*xShift).^2+(v(:,2)-ones(nLen,1)*yShift).^2).^0.5;
%% SPC
elseif strcmp(method,'SPC')
    v = trackStackFlow(cat(3,refimg,beadimg),beads,...
        corLength,corLength,'maxSpd',maxFlowSpeed,...
        'mode','accurate');
    err = ((v(:,1)-uBeads(:,1)).^2+(v(:,2)-uBeads(:,2)).^2).^0.5;
else
    error('the method argument should be one of CDWS, CCWS, and SPC');
end

%% Uncertainty estimation (rms)
end
% errCDWS = ((vCDWS(:,1)-uBeads(:,1)).^2+(vCDWS(:,2)-uBeads(:,2)).^2).^0.5;
% errCCWS = ((vCCWS(:,1)-uBeads(:,1)).^2+(vCCWS(:,2)-uBeads(:,2)).^2).^0.5;
% errSPC = ((vSPC(:,1)-uBeads(:,1)).^2+(vSPC(:,2)-uBeads(:,2)).^2).^0.5;
% nLocBeads = length(beads);
% meanErrCDWS = nanmean(errCDWS);
% meanErrCCWS = nanmean(errCCWS);
% meanErrSPC = nanmean(errSPC);