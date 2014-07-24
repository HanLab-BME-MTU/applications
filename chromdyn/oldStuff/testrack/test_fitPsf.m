function results = test_fitPsf
%function to test the influence of different movie conditions on the psf-fitting

%define parameters
positions = [[32.5;32.5;9.5]];% [32;32;9.5] [32.5;32.5;9]

SNRs = [1000000];

backgrounds = [0];

psfs = [2];

patchSizeFactors = [2 3];

backgroundEstimates = [4];

FT_SIGMA = [1.5289 1.1756];

numOfRepeats = 1;

results.mean = [];
results.std = [];

maxval = length(psfs)*size(positions,2)*length(backgrounds)*length(backgroundEstimates)*...
    length(patchSizeFactors)*numOfRepeats*length(SNRs);
waitbarHandle=mywaitbar(0,[],maxval,'testing');
counter = 0;

for PSF = psfs
    options.psf = PSF;
    for POS = positions
        slist(1).sp(1).cord =POS';
        options.movie = 1;
        options.background = -1;
        options.save = 0;
        [synthMovie, dummy, dataProperties, dummy, sMBG] = ...
            generateSynthMovie(slist,[],[],options);
        for BGE = backgroundEstimates
        for SNR = SNRs
        for BG = backgrounds
            for PATCH = patchSizeFactors
                
                    
                        for i = 1:numOfRepeats
                            options.movie = 2;
                            options.background = BG;
                            [testImg,dummy,dummy,dummy,realImgBG(i)] = generateSynthMovie([],[],SNR,options,synthMovie);
                            
                            %select patch. make sure it will be an odd size
                            patchCenter = [32.5,32.5,9.5];
                            delta = ceil(PATCH * FT_SIGMA);
                            %don't correct for possible out of bounds
                            %minBounds = max(floor(patchCenter-delta([1 1 2])),[0 0 0]);
                            %maxBounds = min(ceil(patchCenter+delta([1 1 2])),[64 64 18]);
                            minBounds = round(patchCenter)-delta([1 1 2]);
                            maxBounds = round(patchCenter)+delta([1 1 2]);
                            
                            
                            imgPatch = squeeze(testImg(minBounds(1):maxBounds(1),minBounds(2):maxBounds(2),minBounds(3):maxBounds(3),1,1));
                            
                            switch BGE
                                case 1
                                    
                                    imm = mean(testImg(:));
                                    estBackground(i) = imm;
                                case 2
                                    minImg = min(testImg(:));
                                    for j = 1:18,nse(j)=imNoiseEstim(testImg(:,:,j));end
                                    estBackground(i) = 4*mean(nse)+minImg;
                                case 3
                                    estBackground(i)=realImgBG(i);
                                case 4
                                    imm = median(testImg(:));
                                    estBackground(i) = imm;
                                case 5
                                    imm = mean(testImg(:));
                                    estBackground(i) = 1.1*imm;
                            end
                            
                            %prepare initial parameters (dx, dy, dz, amp (not used in the fit), sigXY, sigZ)
                            parms1 = [0.5,0.5,0.5, 0, FT_SIGMA(1), FT_SIGMA(2)];
                            
                            [dummy,fittedParms(i,:), parmGradient]=fitGauss2PSF(imgPatch,parms1,estBackground(i));
                            
                            counter = counter+1;
                            mywaitbar(counter/maxval,waitbarHandle,maxval);
                        end
                        fittedParms(:,1:3)=fittedParms(:,1:3)/0.5;
                        fittedParms(:,4)=1;
                        fittedParms(:,5)=fittedParms(:,5)/FT_SIGMA(1);
                        fittedParms(:,6)=fittedParms(:,6)/FT_SIGMA(2);
                        mP = mean(fittedParms,1);
                        stdP = std(fittedParms,1)./mP;
                        relBGdelta = mean(estBackground./realImgBG);
                        relBGstd = std(estBackground./realImgBG)/relBGdelta;
                        results.mean = [results.mean;mP(1),mP(2),mP(3),mP(5),mP(6),relBGdelta,SNR,PATCH,BG,BGE];
                        results.std = [results.std;stdP(1),stdP(2),stdP(3),stdP(5),stdP(6),relBGstd,SNR,PATCH,BG,BGE];
                        
                    end %for SNR = SNRs
                end %for BGE = backgroundEstimates
            end %for PATCH = patchSizeFactors
        end %for BG = backgrounds
    end %for POS = positions
end %for PSF = psfs

close(waitbarHandle)