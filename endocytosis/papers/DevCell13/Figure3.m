%=========================================================================================
% Figure 3
%=========================================================================================
% This script generates the panels for Figure 3 from Aguet et al., Dev. Cell, 2013.

% set to true to print figures
printFigures = false;
f3path = '/Users/aguet/Documents/Papers/DevCell13/Figure 3 - Visitors/';

% root = '/home/fa48/lccb-endocytosis/';
root = '/home/fa48/Documents/MATLAB/Data/';
% dpath = [root 'Marcel_fs003/110301_BSC1_rat_monkey/BSC1_LCa_ratox/'];
dpath = [root '110301_BSC1_rat_monkey/BSC1_LCa_ratox/'];
dataOX = loadConditionData(dpath, {''}, {'EGFP'});

fset = loadFigureSettings('print');
%%
%==================================================================================
% Panel A: Intensity cohorts, for all & above threshold tracks
%==================================================================================
% Initially shown for a single cell (median intensity): Cell5_1s or Cell2!_1s
cohortBounds = [10 20 30 40 60 80 100 120];
xa = (cohortBounds(1:end-1)+cohortBounds(2:end))/2;
rPos = [-2 55 38 60];

opts = {'ShowBackground', false, 'DisplayMode', 'print', 'YTick', 0:25:175,...
    'lftDataName', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat', 'Rescale', false,...
    'FillMode', 'SEM', 'CohortBounds_s', cohortBounds, 'FrontLayer', false,...
    'ExcludeVisitors', false};

[cThr, rThr] = plotIntensityCohorts(dataOX, opts{:}, 'MaxIntensityThreshold', 100,...
    'DisplayMode', 'print', 'YTick', 0:25:175, 'YLim', [0 125]);
set(gca, 'YLim', [0 125], 'XTick', xa, 'XLim', [-10 120]);
rectangle('Position', rPos, 'LineWidth', 0.75);
if printFigures
    print('-depsc2', '-loose', ['IntensityCohorts_dataOX_allCells_CCPs.eps']);
end

% aw = rPos(1)+rPos(3); %130 x-range
% ah = rPos(2)+rPos(4); %125 y-range
% ahc = rPos(4)/125*3.5;
% awc = rPos(3)/130*6;
% 
% % copy figure
% ha = gca;
% h1 = gcf;
% h2 = figure(fset.fOpts{:}, 'Position', get(h1, 'Position'));
% objects = allchild(h1);
% copyobj(get(h1,'children'),h2);
% %%
% axis([rPos(1) aw rPos(2) ah]);
% aw = awc/ahc*3.5;
% set(gca, 'Position', [0.7 2 aw 3.5], 'TickLength', fset.TickLength*6/aw);
% xt = get(ha, 'XTick')-1;
% set(gca, 'YTick', 50:10:130, 'XTick', []);
% xlabel([]);
% ylabel([]);
% if printFigures
%     print('-depsc2', '-loose', ['IntensityCohorts_dataOX_' getCellDir(dataOX(k)) '_CCPs_inset.eps']);
% end

%%
%==================================================================================
% Panel B: visitor track examples
%==================================================================================
i = 2; % Cell2!_1s
lftDataOX = getLifetimeData(dataOX(i), 'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat',...
    'Overwrite', false, 'ExcludeVisitors', false, 'Cutoff_f', 5, 'ReturnValidOnly', true);
vidx = getVisitorIndex(lftDataOX);
tracks = loadTracks(dataOX(i), 'Category', 'all', 'Mask', true, 'Cutoff_f', 5,...
    'FileName', 'ProcessedTracks_gap=2_radius=5,10_a=2.5.mat');
% cmeDataViewer(dataOX(i), 'LoadFrames', false, 'Cutoff_f', 5);


fset = loadFigureSettings('print');
candIdx = [11845 12393 14164];
candTrue = [11826 12777 14092];
XTick = 0:5:120;

dt = 1;
for k = 1:numel(candIdx)
    ah = 1.5;
    aw = 3;
    itrack = tracks(candIdx(k));

    wref = fset.axPos(3); % reference
    width = itrack.t(end)-itrack.t(1) + 7*dt + 7*dt;
    xscale = width/30;
    aw = aw*xscale;
    if aw>2
        tickLength = fset.TickLength*wref/aw;
    else
        tickLength = fset.TickLength*wref/2;
    end
    XLim = [-7*dt itrack.t(end)-itrack.t(1)+7*dt];
    
    figure(fset.fOpts{:}, 'Position', [6 6 4 2.5]);
    ha = axes(fset.axOpts{:}, 'Position', [0.75 0.5 aw ah], 'TickLength', tickLength);
    hold on;
    plotTrack(dataOX(i), itrack, 1, 'Handle', ha, 'DisplayMode', 'print',...
        'MarkerSizes', [6 1.5 0.5], 'LineWidth', 0.5);
    set(ha, 'XLim', XLim, 'XTick', XTick, 'YLim', [-50 250], 'YTick', 0:50:300);
    hold on;
    
    itrackRef = tracks(candTrue(k));
    plotTrack(dataOX(i), itrackRef, 1, 'Handle', ha, 'DisplayMode', 'print',...
        'MarkerSizes', [6 1.5 0.5], 'LineWidth', 0.5, 'Background', 'off', 'Hues', [0.6 0 0]);
    if printFigures
        print('-depsc2', '-loose', [f3path 'vTrack_' num2str(candIdx(k), '%.4d') '.eps']);
    end
    
    [stack, xa, ya] = getTrackStack(dataOX(i), itrack,...
        'WindowWidth', 6, 'Reference', 'frame');
    plotTrackMontage(itrack, stack, xa, ya,'DynamicRange', 231+[0 250],...
        'ShowMarkers', false, 'ShowDetection', false, 'FramesPerRow', 20, 'Width', 600);
    if printFigures
        print('-depsc2', '-loose', [f3path 'vMontage_' num2str(candIdx(k), '%.4d') '.eps']);
    end
end

%%
%==================================================================================
% Panel C: lifetime distribution after visitor exclusion
%==================================================================================
lftResOX = runLifetimeAnalysis(dataOX, 'LifetimeData', 'LifetimeData_gap=2_radius=5,10_a=2.5.mat',...
    'MaxIntensityThreshold', 100, 'Cutoff_f', 5, 'Display', 'off');
plotLifetimes(lftResOX, 'DisplayMode', 'print', 'PlotAll', true);
if printFigures
    print('-depsc2', '-loose', [f3path 'lftDist_FinalThresold.eps']);
end
