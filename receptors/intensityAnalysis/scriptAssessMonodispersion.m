%Amplitude percentiles from all frames combined

ampAll10_1 = vertcat(movieInfo10_1.amp); ampAll10_1 = ampAll10_1(:,1);
ampAll10_2 = vertcat(movieInfo10_2.amp); ampAll10_2 = ampAll10_2(:,1);
ampAll10_3 = vertcat(movieInfo10_3.amp); ampAll10_3 = ampAll10_3(:,1);

ampAll20_1 = vertcat(movieInfo20_1.amp); ampAll20_1 = ampAll20_1(:,1);
ampAll20_2 = vertcat(movieInfo20_2.amp); ampAll20_2 = ampAll20_2(:,1);
ampAll20_3 = vertcat(movieInfo20_3.amp); ampAll20_3 = ampAll20_3(:,1);

ampAll40_1 = vertcat(movieInfo40_1.amp); ampAll40_1 = ampAll40_1(:,1);
ampAll40_2 = vertcat(movieInfo40_2.amp); ampAll40_2 = ampAll40_2(:,1);
ampAll40_3 = vertcat(movieInfo40_3.amp); ampAll40_3 = ampAll40_3(:,1);

prctile00_10 = [prctile(ampAll10_1,0) prctile(ampAll10_2,0) prctile(ampAll10_3,0)];
prctile00_20 = [prctile(ampAll20_1,0) prctile(ampAll20_2,0) prctile(ampAll20_3,0)];
prctile00_40 = [prctile(ampAll40_1,0) prctile(ampAll40_2,0) prctile(ampAll40_3,0)];

prctile01_10 = [prctile(ampAll10_1,1) prctile(ampAll10_2,1) prctile(ampAll10_3,1)];
prctile01_20 = [prctile(ampAll20_1,1) prctile(ampAll20_2,1) prctile(ampAll20_3,1)];
prctile01_40 = [prctile(ampAll40_1,1) prctile(ampAll40_2,1) prctile(ampAll40_3,1)];

prctile99_10 = [prctile(ampAll10_1,99) prctile(ampAll10_2,99) prctile(ampAll10_3,99)];
prctile99_20 = [prctile(ampAll20_1,99) prctile(ampAll20_2,99) prctile(ampAll20_3,99)];
prctile99_40 = [prctile(ampAll40_1,99) prctile(ampAll40_2,99) prctile(ampAll40_3,99)];

prctile00 = ...
    [mean(prctile00_10) mean(prctile00_20) mean(prctile00_40); ...
    std(prctile00_10) std(prctile00_20) std(prctile00_40)]';

prctile01 = ...
    [mean(prctile01_10) mean(prctile01_20) mean(prctile01_40); ...
    std(prctile01_10) std(prctile01_20) std(prctile01_40)]';

prctile99 = ...
    [mean(prctile99_10) mean(prctile99_20) mean(prctile99_40); ...
    std(prctile99_10) std(prctile99_20) std(prctile99_40)]';

%Number of particles per frame

for i=1:500, numPart10_1(i) = size(movieInfo10_1(i).xCoord,1); end
for i=1:500, numPart10_2(i) = size(movieInfo10_2(i).xCoord,1); end
for i=1:500, numPart10_3(i) = size(movieInfo10_3(i).xCoord,1); end

for i=1:500, numPart20_1(i) = size(movieInfo20_1(i).xCoord,1); end
for i=1:500, numPart20_2(i) = size(movieInfo20_2(i).xCoord,1); end
for i=1:500, numPart20_3(i) = size(movieInfo20_3(i).xCoord,1); end

for i=1:500, numPart40_1(i) = size(movieInfo40_1(i).xCoord,1); end
for i=1:500, numPart40_2(i) = size(movieInfo40_2(i).xCoord,1); end
for i=1:500, numPart40_3(i) = size(movieInfo40_3(i).xCoord,1); end

numPartFirst5frames10 = [mean(numPart10_1(1:5)) mean(numPart10_2(1:5)) mean(numPart10_3(1:5))];
numPartFirst5frames20 = [mean(numPart20_1(1:5)) mean(numPart20_2(1:5)) mean(numPart20_3(1:5))];
numPartFirst5frames40 = [mean(numPart40_1(1:5)) mean(numPart40_2(1:5)) mean(numPart40_3(1:5))];

numPartFirst5frames = ...
    [mean(numPartFirst5frames10) mean(numPartFirst5frames20) mean(numPartFirst5frames40); ...
    std(numPartFirst5frames10) std(numPartFirst5frames20) std(numPartFirst5frames40)]';

%Integrated intensity per particle

for i=1:500, movieInfo10_1(i).integratedIntensity = movieInfo10_1(i).amp(:,1)*2*pi*psfSigma^2; end
for i=1:500, movieInfo10_2(i).integratedIntensity = movieInfo10_2(i).amp(:,1)*2*pi*psfSigma^2; end
for i=1:500, movieInfo10_3(i).integratedIntensity = movieInfo10_3(i).amp(:,1)*2*pi*psfSigma^2; end

for i=1:500, movieInfo20_1(i).integratedIntensity = movieInfo20_1(i).amp(:,1)*2*pi*psfSigma^2; end
for i=1:500, movieInfo20_2(i).integratedIntensity = movieInfo20_2(i).amp(:,1)*2*pi*psfSigma^2; end
for i=1:500, movieInfo20_3(i).integratedIntensity = movieInfo20_3(i).amp(:,1)*2*pi*psfSigma^2; end

for i=1:500, movieInfo40_1(i).integratedIntensity = movieInfo40_1(i).amp(:,1)*2*pi*psfSigma^2; end
for i=1:500, movieInfo40_2(i).integratedIntensity = movieInfo40_2(i).amp(:,1)*2*pi*psfSigma^2; end
for i=1:500, movieInfo40_3(i).integratedIntensity = movieInfo40_3(i).amp(:,1)*2*pi*psfSigma^2; end

%Sum of integrated intensity 

for i=1:500, sumIntegIntPerFrame10_1(i) = sum(movieInfo10_1(i).integratedIntensity); end
for i=1:500, sumIntegIntPerFrame10_2(i) = sum(movieInfo10_2(i).integratedIntensity); end
for i=1:500, sumIntegIntPerFrame10_3(i) = sum(movieInfo10_3(i).integratedIntensity); end

for i=1:500, sumIntegIntPerFrame20_1(i) = sum(movieInfo20_1(i).integratedIntensity); end
for i=1:500, sumIntegIntPerFrame20_2(i) = sum(movieInfo20_2(i).integratedIntensity); end
for i=1:500, sumIntegIntPerFrame20_3(i) = sum(movieInfo20_3(i).integratedIntensity); end

for i=1:500, sumIntegIntPerFrame40_1(i) = sum(movieInfo40_1(i).integratedIntensity); end
for i=1:500, sumIntegIntPerFrame40_2(i) = sum(movieInfo40_2(i).integratedIntensity); end
for i=1:500, sumIntegIntPerFrame40_3(i) = sum(movieInfo40_3(i).integratedIntensity); end

%Sum of integrated intensity divided by number of particles in first frame

sumIntegIntPerFrameDividedByNumPartFrame1_10_1 = sumIntegIntPerFrame10_1/numPart10_1(1);
sumIntegIntPerFrameDividedByNumPartFrame1_10_2 = sumIntegIntPerFrame10_2/numPart10_2(1);
sumIntegIntPerFrameDividedByNumPartFrame1_10_3 = sumIntegIntPerFrame10_3/numPart10_3(1);

sumIntegIntPerFrameDividedByNumPartFrame1_20_1 = sumIntegIntPerFrame20_1/numPart20_1(1);
sumIntegIntPerFrameDividedByNumPartFrame1_20_2 = sumIntegIntPerFrame20_2/numPart20_2(1);
sumIntegIntPerFrameDividedByNumPartFrame1_20_3 = sumIntegIntPerFrame20_3/numPart20_3(1);

sumIntegIntPerFrameDividedByNumPartFrame1_40_1 = sumIntegIntPerFrame40_1/numPart40_1(1);
sumIntegIntPerFrameDividedByNumPartFrame1_40_2 = sumIntegIntPerFrame40_2/numPart40_2(1);
sumIntegIntPerFrameDividedByNumPartFrame1_40_3 = sumIntegIntPerFrame40_3/numPart40_3(1);


%Sum of integrated intensity divided by number of particles in first frame
%divided by background std

sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_10_1 = sumIntegIntPerFrameDividedByNumPartFrame1_10_1/bkgStd10;
sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_10_2 = sumIntegIntPerFrameDividedByNumPartFrame1_10_2/bkgStd10;
sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_10_3 = sumIntegIntPerFrameDividedByNumPartFrame1_10_3/bkgStd10;

sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_20_1 = sumIntegIntPerFrameDividedByNumPartFrame1_20_1/bkgStd10;
sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_20_2 = sumIntegIntPerFrameDividedByNumPartFrame1_20_2/bkgStd10;
sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_20_3 = sumIntegIntPerFrameDividedByNumPartFrame1_20_3/bkgStd10;

sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_40_1 = sumIntegIntPerFrameDividedByNumPartFrame1_40_1/bkgStd10;
sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_40_2 = sumIntegIntPerFrameDividedByNumPartFrame1_40_2/bkgStd10;
sumIntegIntPerFrameDivByNumPartFrame1DivByBkgStd_40_3 = sumIntegIntPerFrameDividedByNumPartFrame1_40_3/bkgStd10;



