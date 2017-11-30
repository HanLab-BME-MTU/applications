% Assaf June. 2017
function [] = control_focus()

close all;

% addpath(genpath('/home2/azaritsky/code/common'));
% addpath(genpath('/home2/azaritsky/code/extern'));
% 
% % addpath(genpath('/home2/azaritsky/code/applications/monolayer/utils'));
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/algs'));
% addpath(genpath('/home2/azaritsky/code/applications/monolayer/timeLapseAnalysis'));
addpath(genpath('/home2/azaritsky/code/applications/2dActionRecognition'));
addpath(genpath('/home2/azaritsky/code/applications/monolayer/utils/'));

dataDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/raw/Controls/Focus_change/';
analysisDir = '/project/bioinformatics/Danuser_lab/liveCellHistology/analysis/Controls/Focus_change/';
outdir = [analysisDir filesep 'out' filesep];

if ~exist(outdir,'dir')
    mkdir(outdir);
end

always = false;
experiments = {{'190613_m481_unfocus001.nd2','190613_m481_focus.nd2','190613_m481_unfocus002.nd'}};
ntasks = 3;
nexperiments = length(experiments);

N = ntasks * nexperiments;

lbpMapping = getmapping(8,'riu2');

nScale = 4;

lbpAll = cell(1,nScale);
for iScale = 1 : nScale
    lbpAll{iScale}.lbp = nan(10,N);
end

for iexp = 1 : nexperiments
    curExp = experiments{iexp};
    
    unfocus1 = curExp{1};
    focus = curExp{2};
    unfocus2 = curExp{3};
        
    [~, nameUnfocus1, ~] = fileparts(unfocus1);
    [~, nameFocus, ~] = fileparts(focus);
    [~, nameUnfocus2, ~] = fileparts(unfocus2);
  
    %% Create movie data
    if ~exist([analysisDir nameUnfocus1],'dir') && ~always
        try
            MDUnfocus1 = MovieData(fullfile([dataDir nameUnfocus1 '.nd2']),'outputDirectory', fullfile([analysisDir nameUnfocus1 filesep nameUnfocus1]));
            MDFocus = MovieData(fullfile([dataDir nameFocus '.nd2']),'outputDirectory', fullfile([analysisDir nameFocus filesep nameFocus]));
            MDUnfocus2 = MovieData(fullfile([dataDir nameUnfocus2 '.nd2']),'outputDirectory', fullfile([analysisDir nameUnfocus2 filesep nameUnfocus2]));
        catch ee
            warning(['file ' [dataDir name '.nd2'] 'can not be opend by MovieData']);
            continue;
        end
    end
    
    %%
    for itask = 1 : ntasks
        MDUnfocus1Fname = [analysisDir nameUnfocus1 filesep nameUnfocus1 '_s' num2str(itask) filesep nameUnfocus1 '_s' num2str(itask) '.mat'];
        MDFocusFname = [analysisDir nameFocus filesep nameFocus '_s' num2str(itask) filesep nameFocus '_s' num2str(itask) '.mat'];
        MDUnfocus2Fname = [analysisDir nameUnfocus2 filesep nameUnfocus2 '_s' num2str(itask) filesep nameUnfocus2 '_s' num2str(itask) '.mat'];        
        
        MDUnfocus1 = MovieData.load(MDUnfocus1Fname);
        MDFocus = MovieData.load(MDFocusFname);
        MDUnfocus2 = MovieData.load(MDUnfocus2Fname);
        
        IUnfocus1 = MDUnfocus1.getChannel(1).loadImage(1);
        IFocus = MDFocus.getChannel(1).loadImage(1);
        IUnfocus2 = MDUnfocus2.getChannel(1).loadImage(1);        
        
        imwrite(IUnfocus1,[outdir filesep nameUnfocus1 '_s' num2str(itask) '.tif'],'tif');
        imwrite(IFocus,[outdir filesep nameFocus '_s' num2str(itask) '.tif'],'tif');
        imwrite(IUnfocus2,[outdir filesep nameUnfocus2 '_s' num2str(itask) '.tif'],'tif');
        
        figure; imagesc(IUnfocus1); colormap(gray); saveas(gcf,[outdir filesep nameUnfocus1 '_s' num2str(itask) '.jpg']);
        figure; imagesc(IFocus); colormap(gray); saveas(gcf,[outdir filesep nameFocus '_s' num2str(itask) '.jpg']);
        figure; imagesc(IUnfocus2); colormap(gray); saveas(gcf,[outdir filesep nameUnfocus2 '_s' num2str(itask) '.jpg']);
        
        Ireg = zeros(size(IUnfocus1,1),size(IUnfocus1,2),3);
        Ireg(:,:,1) = (IUnfocus1 - mean(IUnfocus1(:)))./std(double(IUnfocus1(:)));
        Ireg(:,:,2) = (IFocus - mean(IFocus(:)))./std(double(IFocus(:)));
        Ireg(:,:,3) = (IUnfocus2 - mean(IUnfocus2(:)))./std(double(IUnfocus2(:)));
        
        imwrite(Ireg,[outdir filesep nameUnfocus1(1:11) '_s' num2str(itask) '.tif'],'tif');
        
        figure; imagesc(Ireg); saveas(gcf,[outdir filesep nameUnfocus1(1:11) '_s' num2str(itask) '.jpg']);
        close all;
        %         [optimizer, metric] = imregconfig('monomodal');
        %         optimizer.MaximumStepLength = 3;
        %         I10x_zoom_reg = imregister(I10x_zoom,I20x10x,'translation',optimizer, metric);
        
        lbpTFname = [outdir filesep nameUnfocus1 '_s' num2str(itask) '.mat'];
        
        IUnfocus1_pyramidLBP = getPyramidLbp(IUnfocus1,nScale,lbpMapping);
        IFocus1_pyramidLBP = getPyramidLbp(IFocus,nScale,lbpMapping);
        IUnfocus2_pyramidLBP = getPyramidLbp(IUnfocus2,nScale,lbpMapping);
        
        pyramidLBP.unfocus1 = IUnfocus1_pyramidLBP;
        pyramidLBP.focus = IFocus1_pyramidLBP;
        pyramidLBP.unfocus2 = IUnfocus2_pyramidLBP;
        
        for iScale = 1 : nScale            
            lbpAll{iScale}.lbp(:,(itask-1)*ntasks+1) = pyramidLBP.unfocus1{iScale}.lbpDesc;
            lbpAll{iScale}.lbp(:,(itask-1)*ntasks+2) = pyramidLBP.focus{iScale}.lbpDesc;
            lbpAll{iScale}.lbp(:,(itask-1)*ntasks+3) = pyramidLBP.unfocus2{iScale}.lbpDesc;
        end
        
        save(lbpTFname,'pyramidLBP','lbpAll');
    end        
end

lbpSimilarity = cell(1,4);
for iScale = 1 : nScale
    features = lbpAll{iScale}.lbp';
    
    lbpSimilarity{iScale} = squareform(pdist(features,'cityblock'));        
    
    h = figure; hold on; imagesc(lbpSimilarity{iScale}); colorbar; hold off; saveas(h,[outdir filesep 'sim_' num2str(iScale) '.tif']);
    
    [nObs mFeats] = size(features);
    meanFeats = mean(features);
    stdFeats = std(features);
    curFeatsNorm = (features - repmat(meanFeats,[nObs 1])) ./ repmat(stdFeats,[nObs 1]); % zscore(A)
    
    [coeff,score,latent] = pca(curFeatsNorm);
    
    accVariance = cumsum(latent)./sum(latent);
    
    plotPCA(score,[outdir filesep 'pca_' num2str(iScale) '.tif']);
    
    close all;
end
end

function [] = plotPCA(score,outFname)
marker = {'s','o','d'};
colors = {'r','g','b'};

h = figure;
hold on;
for itask = 1 : 3
    for ifocus = 1 : 3
        plot(score((itask-1)*3+ifocus,1),score((itask-1)*3+ifocus,2),sprintf('%s%s',marker{ifocus},colors{itask}),...
            'MarkerSize',8,'MarkerFaceColor',colors{itask});
    end
end
hold off;
saveas(h,outFname);
end

function pyramidLBP = getPyramidLbp(I,nScale,lbpMapping)
pyramidLBP = cell(1,nScale);
for iScale = 1 : nScale
    curScale = 2^(-iScale+1);
    curI = imresize(I,curScale);
    
    IlbpTmp = lbp(curI,1,8,lbpMapping,'');
    Ilbp = nan(size(curI));
    Ilbp(2:end-1,2:end-1) = IlbpTmp;
    pyramidLBP{iScale}.Ilbp = Ilbp;
    
    lbpDesc = hist(Ilbp(:),0:9);
    pyramidLBP{iScale}.lbpDesc = lbpDesc ./ sum(lbpDesc);
end
end
	
	
