if(~exist('I','var'))
    [~,hostname] = system('hostname');
    hostname = strtrim(hostname);

    switch(hostname)
        case 'FSM2UA220003Q'
            % HP Z820 Goldman workstation
            cd 'Z:\Takeshi\N-SIM\040715';
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.nd2');
            I = MD.channels_(1).loadImage(1,11);
        case 'mkitti-jaqaman'
            % Laptop T440s
            load('C:\Users\Mark Kittisopikul\Documents\Data\Lamins\MEFLB1-LACLB12-006_Reconstructed_study\MEFLB1-LACLB12-006_Reconstructed\MEFLB1-LACLB12-006_Reconstructed.mat');
            MD.sanityCheck;
            I = MD.channels_(1).loadImage(1,11);
        otherwise
            % BioHPC
            cd ~/shortcuts/MEFLB1-LACLB12-006_Reconstructed/
            MD = MovieData.load('MEFLB1-LACLB12-006_Reconstructed.mat');
            I = MD.channels_(1).loadImage(1,11);
    end
end
% I = imread('example.tif');
F = OrientationSpaceFilter.constructByRadialOrder(1/2/pi./2,1,8,'none');
R = F*I;
K = 8:-0.1:1;
% t = linspace(1/(2*8+1).^2,1/(2+1).^2,1000);
% K = 1/2./sqrt(t)-1/2;
% rho = zeros(17,length(K));
% for i=1:length(K)
% %     out(:,i) = interpft_extrema(R.getResponseAtOrderFTatPoint(628,323,K(i)));
%       rho(:,i) = R.getResponseAtOrderFTatPoint(628,323,K(i));
% end
% rho = R.getResponseAtOrderFTatPoint(623,383,K);

%% Select point in image
r = 622;
c = 364;
rho = R.getResponseAtOrderFTatPoint(r,c,K);
demo = mat2gray(rho(:,1));

cm = lines;

%% Make figure
figure;
numSamples = 1000;
% findpeaks((interpft(demo,360)-1)*100,(0:359)/2)
plot((0:360-1)/2,(interpft(demo,360)-1)*100,'k');
hold on;
% plot((0:4)/5*180,(interpft(demo,5)-1)*100,'sk');
allSampledMaximaCell = cell(1,numSamples);
sampledTime = zeros(1,numSamples);
for i=1:numSamples
    tic;
    sampled = interpft(demo,i);
    sampledMaxima = (locmax1d(sampled([end 1:end 1]))-2)/i*180;
    sampledMaxima(sampledMaxima < 0) = sampledMaxima(sampledMaxima < 0)+180;
    sampledTime(i) = toc;
    sampledMaxima(:,2) = i;
    allSampledMaximaCell{i} = sampledMaxima;
%     plot(sampledMaxima,repmat(i,length(sampledMaxima),1),'.k');
%     hold on;
end;
allSampledMaxima = vertcat(allSampledMaximaCell{:});
scatter(allSampledMaxima(:,1),allSampledMaxima(:,2),'k.');
ylim([-120 numSamples]);
tic;
calcMaxima = interpft_extrema(demo)/2/pi*180;
calcMaxima = sort(calcMaxima(~isnan(calcMaxima)));
calcTime = toc;
hp = plot([calcMaxima calcMaxima].',repmat([-100 numSamples],length(calcMaxima),1).','--');
for i=1:length(hp)
    set(hp(i),'Color',cm(i,:));
    plot(calcMaxima(i),(interpft1([0 180],demo,calcMaxima(i),'horner')-1)*100,'o','MarkerSize',6,'MarkerEdgeColor',cm(i,:));
end
% hp2 = plot(calcMaxima,(interpft1([0 180],demo,calcMaxima,'horner')-1)*100,'o','MarkerSize',6,'MarkerEdgeColor','r');
grid on;
set(gca,'YTick',0:50:numSamples);
xlim([-1 180]);
xlabel('Orientation (degrees)');
ylabel('Number of Samples');

%% Zoomed in
figure;
% scatter(allSampledMaxima(end-29:end,1),allSampledMaxima(end-29:end,2),'ko');
% hold on;
% hp = plot([calcMaxima calcMaxima].',repmat([numSamples-10 numSamples],length(calcMaxima),1).','--');
allSampledMaxima = vertcat(allSampledMaximaCell{end-10:end});
for i=1:length(calcMaxima)
    subplot(1,3,i);
    scatter(allSampledMaxima(:,1),allSampledMaxima(:,2),'ko');
    hold on;
    plot([calcMaxima(i) calcMaxima(i)].',repmat([numSamples-10 numSamples],1,1).','--','Color',cm(i,:));
    xlim(calcMaxima(i)+[-0.1 0.1]);
    ylim([numSamples-10 numSamples]);
    xlabel('Orientation (degrees)');
ylabel('Number of Samples');
end


%% Show image
figure; imshow(I(r+(0:128)-64,c+(0:128)-64),[]);hold on;
scatter(65,65,'mo')
ylim([33 97]);
xlim([33 97]);
hp = plot(65+bsxfun(@times,[-1 1],sind(calcMaxima)*10).',65+bsxfun(@times,[-1 1],cosd(calcMaxima)*10).','r--');
for i=1:length(hp)
    set(hp(i),'Color',cm(i,:));
%     plot(calcMaxima(i),(interpft1([0 180],demo,calcMaxima(i),'horner')-1)*100,'o','MarkerSize',6,'MarkerEdgeColor',cm(i,:));
end

%% One run Time
figure;
plot(sampledTime);
hold on;
plot([1 numSamples],[calcTime calcTime]);
title('Single Run Time');
xlabel('Number of Samples');
ylabel('Time (s)');

%% Sum absolute error
startNumSampleError = 6;
absError = cellfun(@(x) abs(diff(orientationSpace.diffusion.alignExtrema([x(:,1) calcMaxima],180),1,2)),allSampledMaximaCell(startNumSampleError:numSamples),'UniformOutput',false);
absError = cellfun(@(x) min(180-x,x),absError,'UniformOutput',false);
sumAbsError = cellfun(@(x) sum(x),absError);
figure;
plot(startNumSampleError:numSamples,sumAbsError);
xlabel('Number of Samples');
ylabel('Total Absolute Error (degrees)');
title('Sum Absolute Maxima Error');
grid on;

%% Timeit time calculation
f = @(sampled) (locmax1d(sampled([end 1:end 1]))-2)/i*180;
sampledTimeIt = zeros(1,numSamples);
for i=1:numSamples
    sampledTimeIt(i) = timeit(@() f(interpft(demo,i)));
end;
calcTime = timeit(@() interpft_extrema(demo));

%% Timeit plot
figure;
title('Multiple Run Time using timeit');
xlabel('Number of Samples');
ylabel('Time (s)');
plot(sampledTimeIt);
hold on;
plot([1 numSamples],[calcTime calcTime]);
title('TimeIt Run Time');
xlabel('Number of Samples');
ylabel('Time (s)');
grid on;