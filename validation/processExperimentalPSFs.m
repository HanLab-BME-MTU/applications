%function processExperimentalPSFs


%% ----- Parameters ---- %%

spotDetSig = 1;%Sigma of LoG filter for spot detection
nSigOverBak = 15;%Yeah I know it sounds like overkill but it works.

fitPar = {'X1','X2','X3','A','Sxy','S3','B'};

showPlots = true;

nSigOutlier = 3;%Sigma for detecting outliers in gaussian fits


%% ------------ Load Data ------------- %%

%Parent Directory with separate folders for each PSF image set.
projDir = '/home/he19/files/LCCB/nih/PSFs/Cropped from Cell Images With Beads';
subDirs = dir(projDir);
subDirs = subDirs(3:end);%Get rid of . and .. 

nSub = numel(subDirs);

%Load all the images
for j = 1:nSub    
    imDirs{j} = subDirs(j).name;
    tmp = imDir(imDirs{j});    
    imNames{j} = arrayfun(@(x)(x.name),tmp,'Unif',false);        
    nImages(j) = numel(imNames{j});
    for k = 1:nImages(j)
        
        images{j}{k} = tif3Dread([imDirs{j} filesep imNames{j}{k}]);        
        
    end
    
end

%Just treat them all as independent images
images = [images{:}];
nImagesTot = numel(images);


%% -------- Spot Detection ------------ %%

%Do spot detection on each image

for j = 1:nImagesTot
    
    %Laplacian of gaussian filtering
    logResp = filterLoG(images{j},spotDetSig);
    
    %Local maxima in this response
    locMax{j} = loc_max3Df(logResp);
        
    %Intensity at these local maxima
    locMaxInd{j} = sub2ind(size(images{j}),locMax{j}(:,1),locMax{j}(:,2),locMax{j}(:,3));
    locMaxIntensity{j} = images{j}(locMaxInd{j});
        
    %Estimate image background dist
    [backMean(j),backSTD(j)] = robustMean(double(images{j}(:)));
    
    %Keep local maxima which are much brighter than background   
    locMaxGood{j} = locMaxIntensity{j} > (backMean(j) + nSigOverBak  * backSTD(j));
    locMax{j} = locMax{j}(locMaxGood{j},:);
    locMaxIntensity{j} = locMaxIntensity{j}(locMaxGood{j});
    nSpotsPerIm(j) = nnz(locMaxGood{j});
    
    if 0 %showPlots
        figure
        viewStack(images{j})
        hold on
        plot3(locMax{j}(:,2),locMax{j}(:,1),locMax{j}(:,3),'rx','MarkerSize',10)
    end    
    
end

%% ------------ Gaussian Fitting -------------- %%

%Go through each image and fit a gaussian to the PSF to get sub-pixel
%coordinates of bead

for j = 1:nImagesTot
    
    for k = 1:nSpotsPerIm(j)
    
        %fitPar = {'X1','X2','X3','A','Sxy','Sz','B'};
        initGuess = double([locMax{j}(k,1),locMax{j}(k,2),locMax{j}(k,3),locMaxIntensity{j}(k),1,2,backMean(j)]);        

        [fitGaussPar{j}(k,:),fitGaussSig{j}(k,:),fitGaussCov{j}(k,:,:),fitGaussChi{j}(k),fitGaussDof{j}(k),fitGaussRes{j}{k},fitGaussResGauss{j}{k}] = GaussFitND(double(images{j}),[],fitPar,initGuess);

    end
    
end



%% ------------ Parameter Comparison, Outlier Detection --------- %%


%Combine the parameters from all the fits so we do it per-PSF rather than
%per-image
allFitPar = vertcat(fitGaussPar{:});
allFitSig = vertcat(fitGaussSig{:});
allFitCov = vertcat(fitGaussCov{:});
allFitChi = horzcat(fitGaussChi{:})';
allFitDof = horzcat(fitGaussDof{:})';
allFitRes = [fitGaussRes{:}]';
allFitResGauss =[fitGaussResGauss{:}]';

nSpots = numel(allFitRes);

% TEMP - this doesn't really work because of the aberration (deviation from
% axial symmetry along z) in the experimental PSFs - none of the fits will
% be very good. So since we can't test the residuals for normality, we just
% do outlier detection on the fitted parameters / uncertainties
%Test fits and parameters to exclude shitty ones

shittySpot = false(nSpots,1);


%Test the uncertainties in the localization.
badLoc = detectOutliers(allFitSig(:,1),nSigOutlier);
badLoc = [badLoc detectOutliers(allFitSig(:,2),nSigOutlier)];
badLoc = [badLoc detectOutliers(allFitSig(:,3),nSigOutlier)];
shittySpot(badLoc) = true;

%And the fitted sigma...
badSig = detectOutliers(allFitPar(:,5),nSigOutlier);
badSig = [badSig detectOutliers(allFitPar(:,6),nSigOutlier)];
shittySpot(badSig) = true;

%And amplitude and backround
badAmp = detectOutliers(allFitPar(:,4),nSigOutlier);
badAmp = [badAmp detectOutliers(allFitPar(:,7),nSigOutlier)];
shittySpot(badAmp) = true;



if showPlots
    
    figure
    hist(allFitPar(:,4),20)
    xlabel('Fitted Amplitude')
    ylabel('n')
    
    figure
    hist(allFitPar(:,5),linspace(0,1,20))
    xlabel('Fitted Sigma, X-Y')
    ylabel('n')
    
    figure
    hist(allFitPar(:,6),linspace(0,2,20))
    xlabel('Fitted Sigma, Z')
    ylabel('n')
    
    figure
    hist(allFitPar(:,7),20)
    xlabel('Fitted Background')
    ylabel('n')
    
end

%% ----------- Alignment and Interp ---------- %%


wxy = -7:.25:7;
wz = -7:.25:7;

[xS,yS,zS] = meshgrid(wxy,wxy,wz);%Relative coordinates for interpolation around each PSF

[Mp,Np,Pp] = size(xS);

%Go through each PSF and interpolate image values relative to it's center
for j = 1:nImagesTot
    
    [M,N,P] = size(images{j});
    [X,Y,Z] = meshgrid(1:N,1:M,1:P);
   
    for k = 1:nSpotsPerIm(j)
                
        interpPSFs{j}{k} = interp3(X,Y,Z,double(images{j}), fitGaussPar{j}(k,2) + xS,fitGaussPar{j}(k,1) + yS,fitGaussPar{j}(k,3) + zS);
        %TEMP - WTF??? Why are the returned fitted amplitudes twice the max image
        %intensities??????
        normPSFs{j}{k} = (interpPSFs{j}{k} - fitGaussPar{j}(k,7)) ./  (fitGaussPar{j}(k,4)/2);%Normalize each PSF image to improve averaging.... with an adjustment of a factor of 2 due to some problem with the fitting function... not sure where that came from....??
    end    
    
end



%% ---------- Averaging ------------- %%

%Concatenate since I was using cells to be lazy about initialization and indexing
% allPSFs = [interpPSFs{:}];
% allPSFs = cat(4,allPSFs{:});
allPSFs = [normPSFs{:}];
allPSFs = cat(4,allPSFs{~shittySpot});%Average only those which weren't considered outliers
%allPSFs = cat(4,allPSFs{:});


avgPSF = nanmean(allPSFs,4);
stdPSF = nanstd(allPSFs,[],4);


isoFrac = logspace(log10(.01),log10(.99),10);
nIso = numel(isoFrac);

if showPlots
    
    figure
    subplot(2,2,1)
    %imagesc(squeeze(log(max(avgPSF,[],3))))
    imagesc(squeeze(max(avgPSF,[],3)))
    
    subplot(2,2,2)
    %imagesc(squeeze(log(max(avgPSF,[],2))))
    imagesc(squeeze(max(avgPSF,[],2)))
    
    subplot(2,2,3)
    %imagesc(squeeze(log(max(avgPSF,[],1)))')
    imagesc(squeeze(max(avgPSF,[],1))')
    
    figure
    subplot(2,2,1)
    imagesc(squeeze(max(stdPSF,[],3)))
    
    subplot(2,2,2)
    imagesc(squeeze(max(stdPSF,[],2)))
    
    subplot(2,2,3)
    imagesc(squeeze(max(stdPSF,[],1))')
    
    
    figure
    hold on
    
    cMap = gray(nIso);
    for j = 1:numel(isoFrac)
        fv(j) = isosurface(avgPSF,max(avgPSF(:))*isoFrac(j));
        patch(fv(j),'EdgeColor','none','FaceColor',cMap(j,:),'FaceAlpha',.1);        
    end
    axis equal
    light
    lighting phong
    view(3)
    
end







