function calcMovieFocalAdhesionEnrichment(movieData,varargin)

outDirName = 'adhesion enrichment';%Default name for output directory
outFileName = 'all adhesion enrichment';%File name for exporting to .mat

%% ----------------------- Input ----------------- %%

if nargin < 1 || isempty(movieData)
    [mdFile,mdPath] = uigetfile('*.mat','Please select the MovieData to process:','movieData.mat');
    if mdFile == 0; return, end
    movieData = MovieData.load([mdPath mdFile],0);
    assert(isa(movieData,'MovieData'),'You must select a movieData, not a movieList!');%Because the load is a method of MovieObject....
end

nChanTot = numel(movieData.channels_);

ip = inputParser;

ip.addParamValue('BatchMode',false,@(x)(islogical(x)));
ip.addParamValue('SaveFigures',true,@(x)(islogical(x)));%Whether to save figures to disk
ip.addParamValue('ChannelIndex',1:nChanTot,@(x)(all(ismember(x,1:nChanTot))));
ip.addParamValue('OutputDirectory',[movieData.outputDirectory_ filesep outDirName],@ischar);%Directory for saving figures, spreadsheets etc
ip.addParamValue('MaskChannelIndex',[],@(x)(all(ismember(x,1:nChanTot))));
ip.addParamValue('ColorMap',@jet,@(x)(isa(x,'function_handle') || ischar(x)));%Colormap to use for per-adhesion matrices.
ip.addParamValue('SatPercent',1,@(x)(x >=0 && x < 100));%Percent of pixels to saturate when displaying per-adhesion matrices
ip.addParamValue('Alpha',.01,@(x)(x >=0 && x < 1));%Alpha value for statistical tests
ip.addParamValue('AdhesionColocRadius',0,@(x)(isscalar(x) && x >= 0 && isequal(round(x),x)));%Radius, in pixels, of region around adhesions to include in co-localization analysis

ip.parse(varargin{:});
p = ip.Results;


%% ----------- Init ------------ %%


nChan = numel(p.ChannelIndex);
nFrames = movieData.nFrames_;
tInterval = movieData.timeInterval_;
tData = 0:tInterval:(nFrames-1)*tInterval;
tUnits = 'seconds';
pixSize = movieData.pixelSize_/1e3;
pixUnits = 'um';


if nChan > 2 
    error('This function currently only supports 2 channels!')
end

%Check for segmentation, get the associated process
iSegProc = movieData.getProcessIndex('FocalAdhesionSegmentationProcess',1,~p.BatchMode);
assert(~isempty(iSegProc),'The movie adhesions must be segmented before calcualating adhesion stats!')
segProc = movieData.processes_{iSegProc};

nAdhesions = segProc.maxIndex_;%Total number of adhesions in movie

if isempty(p.MaskChannelIndex)
    iHasMasks = find(segProc.checkChannelOutput(1:nChanTot));
    if nnz(iHasMasks) > 1        
        iSel = listdlg('ListString',arrayfun(@num2str,iHasMasks,'Unif',0),'SelectionMode','single',...
            'Name','Mask Channel Selection','ListSize',[340 314],'PromptString',...
            'Pick channel to use segmented adhesions from:');
        if isempty(iSel); return, end
        p.MaskChannelIndex = iHasMasks(iSel);
    else
        p.MaskChannelIndex = iHasMasks;
    end    
else
    assert(segProc.checkChannelOutput(p.MaskChannelIndex),'No valid masks for selected segmentation process!');
end


%Check the sampling, get the process and output
iSampProc = movieData.getProcessIndex('SegmentationSamplingProcess',1,~p.BatchMode);
assert(~isempty(iSampProc),'The movie segmentation must be sampled before calcualating adhesion enrichment!')
sampProc = movieData.processes_{iSampProc};
assert(all(sampProc.checkChannelOutput(p.ChannelIndex)),'Selected segmentation sampling process does not have valid outputs for selected channel(s)!')
%Make sure that the sampling was done using the mask channel and process
%which was selected
assert(sampProc.funParams_.MaskChannelIndex == p.MaskChannelIndex,'The selected sampling must have been performed using the selected mask channel!')
assert(sampProc.funParams_.SegProcessIndex == iSegProc,'The selected sampling used a different segmentation than the selected segmentation process!')


%Load the sampling
for j = 1:nChan
    allSamples(j) = sampProc.loadChannelOutput(p.ChannelIndex(j));
end
intPropFields = {'MeanIntensity'};
intPropUnits = {'a.u.'};
nIntProp = numel(intPropFields);

mkClrDir(p.OutputDirectory);

varToSave = {'allSamples'};%Array for keeping variables to write to disk. We duplicate the intensity samples for convenience

M = movieData.imSize_(1);
N = movieData.imSize_(2);
P = nFrames;

indChan = p.MaskChannelIndex;%Adhesion-segmented channel is considered independent variable
depChan = p.ChannelIndex(p.ChannelIndex ~= indChan);%Other channel is dependent variable

if p.AdhesionColocRadius > 0
    dStrel = strel('disk',p.AdhesionColocRadius,0);
end

outVars = {'p'};

%% ---------- Image and Mask Loading ------ %%


disp('Loading all images and segmentation...')

allFAMasks = false([M N P]);
allFAColocMasks = false([M N P]);
allIms = cell(nChan,1);%Makes logical indexing easier



for iFrame = 1:nFrames
    
    allFAMasks(:,:,iFrame) = segProc.loadChannelOutput(p.MaskChannelIndex,iFrame); 
    
    if p.AdhesionColocRadius > 0
        allFAColocMasks(:,:,iFrame) = imdilate(allFAMasks(:,:,iFrame),dStrel);
    else
        allFAColocMasks(:,:,iFrame) = allFAMasks(:,:,iFrame);
    end
    
    
    for iChan = 1:nChan
        
        if iFrame == 1
             allIms{iChan} = zeros([M N P],'uint16');
        end
        
       allIms{iChan}(:,:,iFrame) = movieData.channels_(iChan).loadImage(iFrame); 
        
    end
    
end




%% ---------- In-Adhesion Per-Pixel Co-Localization ------ %%

intBins = cell(nChan,1);
binSz = 2;
for j = 1:nChan    
    intBins{j} = min(allIms{j}(:)):binSz:max(allIms{j}(:));
end

intHist2D = hist3([allIms{indChan}(allFAColocMasks(:)) allIms{depChan}(allFAColocMasks(:))],intBins([indChan depChan]));
[lineFit,gof,out] = fit(double(allIms{indChan}(allFAColocMasks(:))),double(allIms{depChan}(allFAColocMasks(:))),'poly1');

fitCoefCI = confint(lineFit,1-p.Alpha);
fitPredCI = predint(lineFit,double(intBins{indChan}),1-p.Alpha);

if ~p.BatchMode
    
    currFig = figure;
    imagesc(intBins{indChan},intBins{depChan},log10(intHist2D'))
    %imagesc(log10(intHist2D))
    axis xy
    saturateImageColormap([],1);
    colormap gray
    hold on    
    y = feval(lineFit,double(intBins{indChan}));
    plot(intBins{indChan},y,'-r','LineWidth',2)
    plot(intBins{indChan},fitPredCI(:,1),'--r','LineWidth',2)
    legend('Fit',['Fit ' num2str( (1-p.Alpha)*100) '% Prediction CI'])
    plot(intBins{indChan},fitPredCI(:,2),'--r','LineWidth',2)
    title({['y = ' num2str(lineFit.p1) '* X + ' num2str(lineFit.p2) ',   Slope ' num2str( (1-p.Alpha)*100) '% CI: ' num2str(fitCoefCI(1,1)) ' to ' num2str(fitCoefCI(2,1))],...
           ['R^2 = ' num2str(gof.rsquare) ',   Mean ' num2str( (1-p.Alpha)*100) '% prediction CI width ' num2str(mean(diff(fitPredCI,1,2)))]})
    xlabel(['Intra-Adhesion Intensity, Channel ' num2str(indChan) ,' a.u.'])
    ylabel(['Intra-Adhesion Intensity, Channel ' num2str(depChan) ,' a.u.'])
    
    if p.SaveFigures
        figName = [p.OutputDirectory filesep 'intra adhesion per pixel colocalization'];
        mfFigureExport(currFig,figName);        
    end
end

outVars = [outVars 'lineFit','gof','out','fitCoefCI','fitPredCI','intBins','intHist2D','indChan','depChan'];


%% ------------- Output --------------- %%

save([p.OutputDirectory filesep outFileName '.mat'],outVars{:})


