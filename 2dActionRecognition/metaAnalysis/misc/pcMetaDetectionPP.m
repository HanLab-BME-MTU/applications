
function [] = pcMetaDetectionPP()

% Usese the motion-estimation cross-correlation scores to calculate
% distribution of the differences between detections and their background.

% Output at analysis/metaAnalysis (jointDistribution)

% Assaf Zaritsky, June 2015

always = true;

addpath(genpath('/home2/azaritsky/code/common/mathfun/psfModels'));
addpath(genpath('/home2/azaritsky/code/common/detectionAlgorithms'));

warning('off','all');

fprintf(sprintf('\nBio format in path %d\n',bfCheckJavaPath()));

analysisDirname = '/project/cellbiology/gdanuser/melanomaModel/Analysis/';
metaDataFname = [analysisDirname 'MetaData/Experiments20150616_AN.mat'];

outFname = [analysisDirname 'metaAnalysis/detectionVsBackgroundScore.mat'];

if exist(outFname,'file') && ~always
    load(outFname);
else
    load(metaDataFname);
    
    detectionMatchScores = [];
    backgroundScore = [];
    for i = 1 : 1 : metaData.tasks.N
        if i > metaData.tasks.N
            return;
        end
        curExp = metaData.tasks.exps(i);
        curTask = metaData.tasks.tasks(i);
        curFname = metaData.experiments.fnames{curExp};
        if curTask <= metaData.experiments.n1{curExp}
            curSource = metaData.experiments.source1{curExp};
        else
            curSource = metaData.experiments.source2{curExp};
        end
        
        %% make sure pcInitiateData was run earlier
        mdFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
            curFname '_s' sprintf('%02d',curTask) filesep...
            curFname '_s' sprintf('%02d',curTask) '.mat'];
        
        if ~exist(mdFname,'file')
            mdFname = [analysisDirname 'Data/' curSource filesep curFname filesep...
                curFname '_s' sprintf('%d',curTask) filesep...
                curFname '_s' sprintf('%d',curTask) '.mat'];
        end
        
        if ~exist(mdFname,'file')
            continue;
        end
        
        MD =  MovieData.load(mdFname);
        for j = 1 : 1000
            detectionPPFname = [MD.outputDirectory_ '/detectCells/detectionsPP/' sprintf('%03d',j) '_detectionStats.mat'];
            
            if ~exist(detectionPPFname,'file')
                continue;
            end
            
            load(detectionPPFname);
            
            detectionMatchScores = [detectionMatchScores detections.LowRes.matchScore];
            backgroundScore = [backgroundScore detections.LowRes.backgroundScore];
        end
    end
    
    diffDetectionBackgroundScore = detectionMatchScores - backgroundScore;
    save(outFname,'detectionMatchScores','backgroundScore','diffDetectionBackgroundScore');
end
    
[jointDistribution] = getScatterQuantification(backgroundScore,detectionMatchScores,0:0.02:1,0:0.02:1);

fontsize = 10;
h = plotJointDistribution(jointDistribution,fontsize);
axisHandle= findobj(h,'type','axes');
set(h,'Color','w');
set(axisHandle,'box','off','XMinorTick','off','TickDir','out','YMinorTick','off','FontSize',fontsize);
set(get(axisHandle,'XLabel'),'FontSize',fontsize); set(get(axisHandle,'YLabel'),'FontSize',fontsize);
export_fig([analysisDirname 'metaAnalysis/jointDistribution.eps']);

figure; hist(diffDetectionBackgroundScore,0:0.02:1);
export_fig([analysisDirname 'metaAnalysis/diffDetectionBackgroundScoreDistribution.eps']);
end

%%
function h = plotJointDistribution(jointDist,fontsize)
h = figure;
imagescnan(jointDist);
hold on;
% caxis(xaxisVals);
haxes = findobj(h,'type','axes');
set(haxes,'XTick',[1,26,51]);
set(haxes,'YTick',[1,26,51]);
set(haxes,'XTickLabel',[0,0.5,1]);
set(haxes,'YTickLabel',[1,0.5,0]);
set(haxes,'FontSize',fontsize);
hold off;
end
