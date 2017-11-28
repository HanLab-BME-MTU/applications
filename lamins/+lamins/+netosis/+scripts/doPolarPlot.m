%% Load Image
% MD = MovieData.load('HL-60DNALB1_10min_012_Reconstructed.mat');
MD = MovieData.load;
CR = CellReader(MD.getReader());

%% Sum channels
lamin = cat(3,CR{2,1,:});
slamin = sum(lamin,3);
I = cat(3,CR{1,1,:});
sI = sum(I,3);

%% Compute images for show
sI_show = imadjust(mat2gray(sI),[0.3 1],[0 1]);
slamin_show = imadjust(mat2gray(slamin),[0.3 1],[0 1]);


%% Acquire mask
% Skip if mask determined in other fashion
askForMask = true;
choices.yes = 'Use Existing Mask';
choices.no = 'Create New Mask';
choices.cancel = 'Cancel';
if(exist('./mask.mat','file'))
    askForMask = questdlg('Use previously created mask?','Mask Exists',choices.yes,choices.no,choices.cancel,choices.yes);
    switch(askForMask)
        case choices.yes
            askForMask = false;
            load('./mask.mat');
        case choices.no
            askForMask = true;
        case choices.cancel
            askForMask = false;
            return;
    end
end

if(askForMask)
    hfig = figure;
    imshowpair(sI_show,slamin_show,'ColorChannels',[2 1 0]);
    h = imellipse;
    wait(h);
    mask = createMask(h);
    close(hfig);
end

%% Get region properties and verify mask
rp = regionprops(mask,'Area','Perimeter','Centroid','BoundingBox','ConvexHull');
hfig = figure; imshowpair(sI_show,slamin_show,'ColorChannels',[2 1 0]);
hold on;
fill(rp.ConvexHull(:,1),rp.ConvexHull(:,2),'w','FaceAlpha',0.3);
saveas(hfig,'mask_overlay.png');
save('mask.mat','mask');

%% Calculate coordinates
[X,Y] = meshgrid(1:size(I,2),1:size(I,1));
XX = X-rp.Centroid(1);
YY = Y-rp.Centroid(2);
[theta,rho] = cart2pol(XX,YY);

offset = (0:16)*2*pi/17;
offset = shiftdim(offset,-1);
theta_offseted = theta - offset;
theta_offseted = wraparoundN(theta_offseted,-pi,pi);
weights = exp(-theta_offseted.^2/2/(2*pi/17).^2).*mask;
normWeights = shiftdim(sum(sum(weights,1),2),2);

%% Process first (DNA) channel
sImasked = sI.*mask;

wsI = shiftdim(weights.*sI,2);
wsI_summed = sum(wsI(:,:),2)./normWeights;

%% Process second (lamin) channel

wslamin = shiftdim(weights.*slamin,2);
wslamin_summed = sum(wslamin(:,:),2)./normWeights;


%% Polar plot
hfig = figure;
imshowpair(sI_show,slamin_show,'ColorChannels',[2 1 0]);
xlim(rp.Centroid(1) + [-300 300]);
ylim(rp.Centroid(2) + [-300 300]);
hold on;
plot(rp.ConvexHull(:,1),rp.ConvexHull(:,2),'y--');

ax2 = polaraxes;
ax2.Color = 'none';
ax2.GridColor = 'w';
ax2.GridAlpha = 1;
ax2.Position = [0.1 0.1 0.8 0.8];
ax2.ThetaColor = 'w';
ax2.RColor = 'w';

hold on;
h = polarplot(mat2gray(interpft(wsI_summed,360)),'g');
h.Parent.ThetaDir = 'clockwise';
h = polarplot(mat2gray(interpft(wslamin_summed,360)),'r');
set(hfig, 'InvertHardCopy', 'off');
saveas(hfig,'polarPlot.fig');
saveas(hfig,'polarPlot.png');

%% Linear plot

hfig = figure; plot(interpft(wsI_summed,360),'g');
xlim([0 360]);
grid on
ax3 = axis;
hold on; plot(interpft(wslamin_summed,360),'r');
hold off
set(gca,'XTick',0:60:360);
set(gca,'XTickLabel',0:60:360);
xlabel('Orientation (Degrees)');
ylabel('Oriented Fluorescence Intensity Summed over Z');
legend({'DNA','Lamin B1'});
saveas(hfig,'linearPlot.fig');
saveas(hfig,'linearPlot.png');





