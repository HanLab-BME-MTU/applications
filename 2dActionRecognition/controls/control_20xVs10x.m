% Assaf June. 2017
function [] = control_20xVs10x()

close all;

% addpath(genpath('/home2/azaritsky/code/common'));
% addpath(genpath('/home2/azaritsky/code/extern'));
% 
% % addpath(genpath('/home2/azaritsky/code/applications/monolayer/utils'));
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/algs'));
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/timeLapseAnalysis'));
% addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));

dataDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/raw/Controls/20x_10x_data/';
analysisDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Controls/20x_10x_data/';
outdir = [analysisDir filesep 'out' filesep];

if ~exist(outdir,'dir')
    mkdir(outdir);
end

always = false;

experiments = {{'170613_m481_10x.nd2','170613_m481_20x.nd2'},{'190613_m481_10x.nd2','190613_m481_20x.nd2'}};

for iexp = 1 : length(experiments)
    curExp = experiments{iexp};
    
    fname10x = curExp{1};
    fname20x = curExp{2};
        
    [~, name10x, ~] = fileparts(fname10x);
    [~, name20x, ~] = fileparts(fname20x);
  
    %% Create movie data
    if ~exist([analysisDir name10x],'dir') && ~always
        try
            MD10x = MovieData(fullfile([dataDir name10x '.nd2']),'outputDirectory', fullfile([analysisDir name10x filesep name10x]));
            MD20x = MovieData(fullfile([dataDir name20x '.nd2']),'outputDirectory', fullfile([analysisDir name20x filesep name20x]));
        catch ee
            warning(['file ' [dataDir name '.nd2'] 'can not be opend by MovieData']);
            continue;
        end
    end
    
    %%
    for itask = 1 : 3
        MD10xFname = [analysisDir name10x filesep name10x '_s' num2str(itask) filesep name10x '_s' num2str(itask) '.mat'];
        MD20xFname = [analysisDir name20x filesep name20x '_s' num2str(itask) filesep name20x '_s' num2str(itask) '.mat'];
        
        MD10x = MovieData.load(MD10xFname);
        MD20x = MovieData.load(MD20xFname);
        
        I10x = MD10x.getChannel(1).loadImage(1);
        [sy,sx] = size(I10x);
        I10x_zoom = I10x(sy/4:(3*sy/4)-1,sx/4:(3*sx/4)-1);
        I20x = MD20x.getChannel(1).loadImage(1);
        I20x10x = imresize(I20x,0.5);
        
        % patch 
        I10x_zoom_shift = I10x_zoom(1:end-34,1:end-22);
        I20x10x_shift = I20x10x(35:end,23:end);
        
        Ireg = zeros(size(I10x_zoom_shift,1),size(I10x_zoom_shift,2),3);
        Ireg(:,:,1) = (I10x_zoom_shift - mean(I10x_zoom_shift(:)))./std(double(I10x_zoom_shift(:)));
        Ireg(:,:,2) = (I20x10x_shift - mean(I20x10x_shift(:)))./std(double(I20x10x_shift(:)));
        
        imwrite(Ireg,[outdir filesep name10x '20x_'  '_s' num2str(itask) '.tif'],'tif');
        
        %         [optimizer, metric] = imregconfig('monomodal');
        %         optimizer.MaximumStepLength = 3;
        %         I10x_zoom_reg = imregister(I10x_zoom,I20x10x,'translation',optimizer, metric);
    end

% pixelSize = 0.325;
% MD =  MovieData.load(mdFname);
end
end

	
	
