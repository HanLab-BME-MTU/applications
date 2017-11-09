function [ corrMatrix ] = getChromaticShift( MD, channels)
%getChromaticShift uses correlation to establish the chromatic shift
%between channels

if(nargin < 2)
    channels = 1:2;
end

import lamins.functions.*;

outputDirectory = [MD.outputDirectory_ filesep 'chromaticShift' filesep];
matFile =  [outputDirectory 'chromaticShift.mat'];

mkdir(outputDirectory);

%% Do calculation
% 1. Calculate correlation matrix
% 2. Calculate means

R = CellReader(CachedReader(MD.getReader()));
if(exist(matFile,'file'))
    load(matFile);
else
    corrMatrix = calcCorrMatrix(R.toCell);
    means = cellfun(@(x) mean(x(:)),R.toCell);
    means = squeeze(means)';
end


%% Plot figures;

h = plotcorrmatrix(corrMatrix,MD.getFilename);
savefig(h,[outputDirectory 'corrMatrix.fig']);

h = figure;
X = squeeze(corrMatrix(channels(1),:,channels(2),:));
imagesc(X);
ylabel(['Channel ' num2str(channels(1))]);
xlabel(['Channel ' num2str(channels(2))]);
th = title(MD.getFilename);
set(th,'interpreter','none');
savefig(h,[outputDirectory 'channelCorr.fig']);

figure;
[~,maxIndex] = max(X);
maxOffset = (1:length(X)) - maxIndex;
bar(maxOffset);
th = title(MD.getFilename);
set(th,'interpreter','none');
xlabel(['Offset: Channel ' num2str(channels(2)) ' - ' num2str(channels(1))]);

savefig(h,[outputDirectory 'offset.fig']);

%% means
h = figure;
nmeans = bsxfun(@minus,means,min(means));
nmeans = bsxfun(@rdivide,nmeans,max(nmeans));
plot(nmeans);
% plot(means);
th = title(MD.getFilename);
set(th,'interpreter','none');
savefig(h,[outputDirectory 'means.fig']);

%% mean diagonal value
h = figure;
X2 = X(10:end-10,10:end-10);
shift = -length(X2):length(X2);
d = zeros(length(shift),1);
idxshift = length(X2) +1;
for i = shift
    d(i+idxshift) = mean(diag(X2,i));
end
plot(shift,d);
[~,maxshiftidx] = max(d);
grid on;
th = title([MD.getFilename ', shift = ' num2str(shift(maxshiftidx))]);
xlabel('Lag: Channel 2 z-position  - Channel 1 z-position');
set(th,'interpreter','none');
savefig(h,[outputDirectory 'diagshift.fig']);

%% xy shift at the z-position
C1 = R(1,1,:).to3D;
C2 = R(2,1,:).to3D;
xc = [];
if(numel(C1) < 1e6)
    xc = pdollar.images.normxcorrn(double(C1),double(C2));
    [~,outmaxi] = max(xc(:));
    [x,y,z] = ind2sub(size(xc),outmaxi);
    delta = [size(xc)+1]/2-[x y z];
    disp(delta);
    
    h = figure;
    th = title([MD.getFilename 'normxcorrn']);
    set(th,'interpreter','none');
    subplot(1,3,1);
    imagesc(squeeze(max(xc,[],3)));
    title('Max Intensity Projection in Z');
    ylabel('y');
    xlabel('x');
    axis ij;
    axis square;
    subplot(1,3,2);
    imagesc(squeeze(max(xc,[],2)));
    title('Max Intensity Projection in X');
    ylabel('y');
    xlabel('z');
    axis ij;
    axis square;
    subplot(1,3,3);
    imagesc(squeeze(max(xc,[],2))');
    title('Max Intensity Projection in Y');
    ylabel('z');
    xlabel('x');
    colormap(jet);
    axis ij;
    axis square;
    savefig(h,[outputDirectory 'xcorrn.fig']);
end

save(matFile,'corrMatrix','X','maxIndex','maxOffset','means','d','xc');

% keyboard;

end

