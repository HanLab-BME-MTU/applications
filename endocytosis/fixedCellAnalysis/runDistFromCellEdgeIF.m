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
ip.addParameter('ChannelNames','', @(x) (ischar(x)||iscell(x)));
ip.parse(data, varargin{:});

%% Adapt data structure for backward compatibility.
for i = 1:numel(data)
    data_old(i).source=data(i).source;
    data_old(i).framePaths=data(i).framePaths;

    for c = 1:length(data(i).channels);
        data_old(i).channels{c} = [data(i).framePaths{c}];
    end
    resDir=[data(i).channels{1} filesep 'Membrane-distance' filesep];
    mkdir(resDir);
    data_old(i).results = [resDir filesep 'Membrane-distance.mat'];
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
analysisFolder=fullfile(fileparts(fileparts(fileparts(data(1).channels{1}))),'analysis');
mkdirRobust(analysisFolder);

fopts = {'Normalized', false, 'Axis', {[0 10 0 100],[0 10 0 1000]},...
    'DisplayMode', 'print','Name',ip.Results.Name, 'Names', {ip.Results.ChannelNames}, ... 
    'Channels', 1:length(data_old(1).channels),'PrintFolder',analysisFolder};
[edgeDistStats,det] = calcDistFromCellEdgeIF(data_old, fopts{:});

for i=1:length(data_old)
    chResults = load(data_old(i).results);

    nCh=length(data_old(i).channels);
    contrastedImg=cell(1,3);
    for c=1:nCh
        chImage = double(imread(data_old(i).channels{c}));
        % build rgb projection (one channel only here)
        dRange1 = [(prctile(chResults.ps(c).c,25)) (prctile(chImage(:), 99.5))];
        contrastedImg{c} = scaleContrast(chImage, dRange1);
    end
    for c=(nCh+1):3
        contrastedImg{c}=zeros(size(contrastedImg{1}));
    end
    rgb = cat(3, contrastedImg{1},contrastedImg{2},contrastedImg{3});

    %% build boundaries
    [ny,nx] = size(chImage);
    B = bwboundaries(chResults.mask);
    B = vertcat(B{:}); % [y x] coordinates
    B(B(:,1)==1 | B(:,1)==ny,:) = [];
    B(B(:,2)==1 | B(:,2)==nx,:) = [];

    %% Overlay
    rgb = uint8(rgb);
    figure('Units', 'Pixels', 'Position', [150 150 nx/2 ny/2], 'PaperPositionMode', 'auto');
    axes('Position', [0 0 1 1]);
    imagesc(rgb); axis image off; colormap(gray(256));
    hold on; plot(B(:,2), B(:,1), 'Color', 0.99*[1 1 1], 'LineWidth', 1);
%    plotScaleBar(5/0.065);
    plotScaleBar(5/0.060,'Label','5 um');
    printPNGEPSFIG(gcf(),analysisFolder,['Cell_' num2str(i) '_mergedChannel']);

end