%runDistFromCellEdgeID(data, varargin) runs CCP detection on a single channel dataset of fixed-cell
% data sets and estimate the distance on the refine mask. The mask is automatically displayed to the 
% for refinement.
% Input:
%    data : structure output of channel data
%
% Options:
%    'Overwrite' : true|{false}
%
% Output:
%     edgeDistStats : structure with fields
%          .dist:           distance to cell edge
%          .A:              amplitude
%          .dfeHists:       distance from edge histogram      
%          .ampHists:       amplitude histogram (as a function of distance)
%          .dfeHistMean_bc: average of bias-corrected dfeHists
%          .dfeHistSD_bc:   s.d. of bias-corrected dfeHists
%                           These are the values plotted
%          .binc:           histogram bin center coordinates
%          .ampHistMean_bc: average of bias-corrected ampHists  
%          .ampHistSD_bc:   s.d. of bias-corrected ampHists
%
%     Note: 'dist' and 'A' are cell arrays of values for each channel
%
%
%     det : structure saved by processFramesIF(), with additional field
%          .ex : coordinates (Nx2) of the cell edge
% Philippe Roudout 12/2017, adapted from FA.

function [edgeDistStats,det]=runDistFromCellEdgeID(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addParameter('Overwrite', false, @islogical);
ip.addParameter('Name', '', @ischar);
ip.parse(data, varargin{:});

%% Adapt data structure for backward compatibility.
nCh=1;
for i = 1:numel(data)
    data_old(i).source=data(i).source;
    data_old(i).framePaths=data(i).framePaths;

    for c = 1:nCh
        data_old(i).channels{c} = [data(i).framePaths{c}];
    end
    mkdir([data(i).channels{c} filesep 'Membrane-distance' filesep]);
    data_old(i).results = [data(i).channels{c} filesep 'Membrane-distance' filesep 'Membrane-distance.mat'];
end


%===============================================================================
% 2) Run detection
%===============================================================================
% reset random number generator to ensure reproducibility
rng('default');
processFramesIF(data_old, 'Overwrite', ip.Results.Overwrite);

%%
%===============================================================================
% 3) Plot 'cell edge distance' histograms
%===============================================================================
fopts = {'Normalized', false, 'Axis', {[0 10 0 100],[0 10 0 1000]},...
    'DisplayMode', 'print','Name',ip.Results.Name, 'Names', {ip.Results.Name},'Channels', [1]};
[edgeDistStats,det] = calcDistFromCellEdgeIF(data_old, fopts{:});


for i=1:length(data_old)
    ch1Image = double(imread(data_old(i).channels{1}));
    ch1Results = load(data_old(i).results);

    %% build rgb projection (one channel only here)
    dRange1 = [(prctile(ch1Results.ps(1).c,25)) (prctile(ch1Image(:), 99.5))];
    tmp1 = scaleContrast(ch1Image, dRange1);
    rgb = cat(3, tmp1, tmp1, tmp1);

    %% build boundaries
    [ny,nx] = size(ch1Image);
    B = bwboundaries(ch1Results.mask);
    B = vertcat(B{:}); % [y x] coordinates
    B(B(:,1)==1 | B(:,1)==ny,:) = [];
    B(B(:,2)==1 | B(:,2)==nx,:) = [];

    %% Overlay
    rgb = uint8(rgb);
    figure('Units', 'Pixels', 'Position', [150 150 nx/2 ny/2], 'PaperPositionMode', 'auto');
    axes('Position', [0 0 1 1]);
    imagesc(rgb); axis image off; colormap(gray(256));
    hold on; plot(B(:,2), B(:,1), 'Color', 0.99*[1 1 1], 'LineWidth', 1);
    plotScaleBar(5/0.065);
end