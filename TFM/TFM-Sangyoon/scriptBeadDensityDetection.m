clear
%% Simulations - with pixel-wise tracking
numBeads = [250 500 750 1000 1500 2000 2500 3000 4000];
nNB = length(numBeads);
pixSize = 108; %nm/px
pixSizeUm = pixSize/1000;
d = 4; %pix (which means ~400 nm)
noiseLevels = 0:0.1:0.2;
nNL = length(noiseLevels);
nExp = 5;

beadDensityDetected = zeros(nNB,nNL,nExp);
w = 2^8; h=w;
fieldArea = w*h; %pix
fieldAreaUm2 = fieldArea*pixSizeUm^2;
beadDensity = numBeads/fieldAreaUm2; % #/um2
B2Bspacing = 1./beadDensity.^2;
%% sigma
bead_d = 40; % nm
pixSize = 108e-9; % m/pix 60x
NA = 1.49; %TIRF
lambda = 665e-9;
M = 1; %1 because I use alreay magnified pixSize.
pixSizeNm=pixSize*1e9;

sigma = getGaussianPSFsigma(NA,M,pixSize,lambda);
%% bead image - single
kon = 100; %s-1
koff = 50; %s-1

nImgs = 100;
padZeros=floor(log10(nImgs))+1;
nDye = 650;
dyeAmp = 1/nDye;

for epm=1:nExp
    p=0;
    ii=0;
    for nB=numBeads
        ii=ii+1;
        jj=0;
        disp(['Currently designed bead density: ' num2str(beadDensity(ii)) '#/\mum^2'])
        for curWN = noiseLevels
            jj=jj+1;
            p=p+1;
            dataPath=['/Volumes/GoogleDrive/My Drive/Documents/Grant/2021 NSF CAREER' ...
                '/Sci/Figs/beadDensitySimulation/SingleFrame/nB' num2str(nB) 'noise' num2str(curWN) ...
                'ex' num2str(epm)];
            if ~exist(dataPath,'dir')
                mkdir(dataPath)
            end
            
            curA = zeros(1,nB);

            % create bead image (single)

            for k=1:nImgs
                for pp=1:nB
                    curA(1,pp) = dyeAmp * sum(rand(1,nDye)<kon/(kon+koff));
                end
                % The simulated fluorescence image produced
                % by this distribution of beads was created by convolution with a Gaussian
                % kernel with σ = 0.21 × λ ,54 where λ is the wavelength of emission NA
                % (here 700 nm) and NA is the numerical aperture of the microscope 
                if k==1 && jj==1
                    [img,bead_x, bead_y, ~, ~] = simGaussianBeads(w, h, sigma, ...
                        'npoints', nB, 'Border', 'truncated','A',curA);
                else
                    img = simGaussianBeads(w, h, sigma, ...
                        'x',bead_x,'y',bead_y,'A', curA, 'Border', 'truncated');
                end
                % camera noise addition
                img = img+curWN*rand(h,w)*max(img(:));
                
                if k==1
                    % save img
                    imwrite(img, [dataPath '.tif'])
                end
                % save img
                imwrite(uint16(img*2^12), [dataPath filesep 'img' num2str(k,['%0.',int2str(padZeros),'d']) '.tif'],'tiff','Compression','none')
            end
    
            if ~exist('hFig','var') || ~isvalid(hFig)
                hFig=figure(1); hold off
            end
            [~,beadDensityDetected(ii,jj,epm)]= ...
                estimateBeadDistance(img,pixSizeNm,sigma,hFig);
            savefig(hFig, [dataPath '.fig'])
            hold off
        end
    end
end
%% save
save('allData.mat')
%% fig
%% plotting
myCs=distinguishable_colors(nNL);
figure, hold on
for jj=1:nNL
    curCol = myCs(jj,:);
    brighterCol = curCol.*1.5;
    brighterCol(brighterCol>1) = 1;
    plot(beadDensity,squeeze(beadDensityDetected(:,jj,:)),'Color',brighterCol);
    h(jj) = plot(beadDensity,mean(beadDensityDetected(:,jj,:),3),'LineWidth',2,'Color',myCs(jj,:));
end
legend(h,'0 %','10%','20%','location','best')

xlabel('Simulated bead density (#/\mum^2)')
ylabel('Detected bead density (#/\mum^2)')
savefig('beadDensityDetected.fig')
%% Do the same thing for LiveSRRF'ed images
% Images are processed via Fiji, LiveSRRF, stored in
% beadDensitySimultation/LiveSRRFed
beadDensityDetectedSRRF = zeros(nNB,nNL,nExp);
mag = 5;
pixSizeNm5x = pixSizeNm/mag;
sigma5x = sigma*mag*0.8;
for epm=1:nExp
    p=0;
    ii=0;
    for nB=numBeads
        ii=ii+1;
        jj=0;
        disp(['Currently designed bead density: ' num2str(beadDensity(ii)) '#/um2'])
        for curWN = noiseLevels
            jj=jj+1;
            p=p+1;
            dataPath=['/Volumes/GoogleDrive/My Drive/Documents/Grant/2021 NSF CAREER' ...
                '/Sci/Figs/beadDensitySimulation/LiveSRRFed/nB' num2str(nB) 'noise' num2str(curWN) ...
                'ex' num2str(epm) ' - liveSRRF (AVG).tif'];
            img = imread(dataPath);
                        
            if ~exist('hFig','var') || ~isvalid(hFig)
                hFig=figure(1); hold off
            end
            [~,beadDensityDetectedSRRF(ii,jj,epm)]= ...
                estimateBeadDistance(img,pixSizeNm5x,sigma5x,hFig);
            savefig(hFig, [dataPath '.fig'])
            hold off
        end
    end
end

%% plotting
myCs=distinguishable_colors(2*nNL);
figure, hold on
for jj=1:nNL
    curCol = myCs(2*jj-1,:);
    brighterCol = curCol.*1.5;
    brighterCol(brighterCol>1) = 1;
    plot(beadDensity,squeeze(beadDensityDetected(:,jj,:)),'Color',brighterCol);
    h(2*jj-1) = plot(beadDensity,mean(beadDensityDetected(:,jj,:),3),'LineWidth',2,'Color',curCol);
    curCol = myCs(2*jj,:);
    brighterCol = curCol.*1.5;
    brighterCol(brighterCol>1) = 1;
    plot(beadDensity,squeeze(beadDensityDetectedSRRF(:,jj,:)),'Color',brighterCol);
    h(2*jj) = plot(beadDensity,mean(beadDensityDetectedSRRF(:,jj,:),3),'LineWidth',2,'Color',curCol);
end
legend(h,'0%-No SR','0%-LiveSRRF','10%-No SR','10%-LiveSRRF','20%-No SR','20%-LiveSRRF','location','best')

xlabel('Simulated bead density (#/\mum^2)')
ylabel('Detected bead density (#/\mum^2)')
savefig('beadDensityDetected.fig')
