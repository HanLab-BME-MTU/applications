xmax=200;
ymax=300;

%% histogram with high density
nPoints = 7000;
sigma = 1.68; % after getGaussianPSFsigma(NA,M,pixSize,lambda);
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
    'npoints', nPoints, 'Border', 'truncated','A',0.3+rand(1,nPoints));
nBeads = length(bead_x);

% randomly placing a NA for 1000 times
nAdh = 5000;
beads = [bead_x bead_y];
d_thres = 5; % pix, assuming f=500 Pa, d=4 pix, distance until u=<0.1

NAx=xmax*rand(nAdh,1);
NAy=ymax*rand(nAdh,1);
[~,dist] = KDTreeClosestPoint(beads, [NAx NAy]);

figure, hist(dist,40)
disp(['density = ' num2str(nBeads/(xmax*ymax/100))])
%% histogram with low density
nPoints = 2000;
[refimg,bead_x, bead_y, ~, Av] = simGaussianBeads(xmax,ymax, sigma, ...
    'npoints', nPoints, 'Border', 'truncated','A',0.3+rand(1,nPoints));
nBeads = length(bead_x);

% randomly placing a NA for 1000 times
nAdh = 5000;
beads = [bead_x bead_y];

NAx=xmax*rand(nAdh,1);
NAy=ymax*rand(nAdh,1);
[~,dist] = KDTreeClosestPoint(beads, [NAx NAy]);

figure, hist(dist,40)
disp(['density = ' num2str(nBeads/(xmax*ymax/100))])
%% plot between bead density vs. percentage of beads experiencing less than .1px of displacement
step=50;
nStep = 20;
d_limit = 8; % f=200Pa, d=8px, u=>0.1 px
pStep=5;
beadDensity = zeros(nStep,1);
percentLostBeads = zeros(nStep,1);
percentLostBeadsErr = zeros(nStep,1);
pLB =  zeros(pStep,1);
for nPoints=step:step:step*nStep

    % randomly placing a NA for 1000 times
    nAdh = 5000;

    for pp=1:pStep
        [~,bead_x, bead_y, ~, ~] = simGaussianBeads(xmax,ymax, sigma, ...
            'npoints', nPoints, 'Border', 'truncated','A',0.3+rand(1,nPoints));
        nBeads = length(bead_x);
        beads = [bead_x bead_y];
        NAx=xmax*rand(nAdh,1);
        NAy=ymax*rand(nAdh,1);
        [~,dist] = KDTreeClosestPoint(beads, [NAx NAy]);
        pLB(pp)=sum(dist>=d_limit)/sum(dist>=0)*100;
    end
    beadDensity(nPoints/step) = nBeads/(xmax*ymax/100); % beads/(10x10 px2)
    percentLostBeads(nPoints/step) = mean(pLB);
    percentLostBeadsErr(nPoints/step) = std(pLB);
end
figure, errorbar(beadDensity,percentLostBeads,percentLostBeadsErr,'k*')
xlabel('Bead density (10-2 #/px2)')
ylabel({'Portion of adhesions whose closest', 'bead is more than 8 px away (%)'})
%% what's the bead density of the real bead image?
beadImgReal = imread('/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/Beads/Cell_11_w2642_t101.TIF');
pstruct = pointSourceDetection(beadImgReal, sigma*0.6, 'alpha', 0.05);
beadDenReal = length(pstruct.x)/(size(beadImgReal,1)*size(beadImgReal,2)/100);
figure, imshow(beadImgReal,[]), hold on, plot(pstruct.x, pstruct.y, 'ro')
%% analyzing real adhesions and forces
forceImgReal = imread('/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/ROI/Colocalization/Region1_gc5/fMap/force100max350.tif');
paxImgReal = imread('/Users/joshua2/Documents/PostdocResearch/Traction Force/Manuscript/Draft1.2/Data/130429_cell11_100f/ROI/Colocalization/Region1_gc5/paxtifs/pax100.tif');
% find the NAs
pstructNA = pointSourceDetection(paxImgReal, 2.1, 'alpha', 0.05);
%find force local maxima in TFM
% fMaxima = locmax2d(forceImgReal,7);
fMaxima = pointSourceDetection(forceImgReal, 3.1, 'alpha', 0.05);
fMaxima =[fMaxima.x' fMaxima.y']; 
NAreal = [pstructNA.x' pstructNA.y'];
[~,distReal] = KDTreeClosestPoint(fMaxima,NAreal);
figure,hist(distReal,20);