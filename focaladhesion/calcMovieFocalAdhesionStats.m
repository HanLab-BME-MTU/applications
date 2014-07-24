function calcMovieFocalAdhesionStats(movieData,varargin)
%CALCMOVIEFOCALADHESIONSTATS calculates focal adhesion stats, makes figures and spreadsheets 
%
% calcMovieFocalAdhesionStats
% calcMovieFocalAdhesionStats(movieData)
%

%Hunter Elliott
%4/2013

%Note - at some point this should be converted into a generic "segmented
%region properties" function, as not much is really focal adhesion
%specifc...


%% ------- Parameters ------- %%

nCols = 64;%Number of colors for colormaps
outDirName = 'adhesion stats';%Default name for output directory
outFileName = 'all adhesion stats';%File name for exporting to .mat

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
assert(~isempty(iSampProc),'The movie segmentation must be sampled before calcualating adhesion stats!')
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


%% ------ Geometric Properties ------ %%

disp('Calculating adhesion geometric properties....')

geoPropFields = {'Area' ,'Eccentricity','Orientation','MajorAxisLength','MinorAxisLength'};%Field names
geoPropUnits = {[pixUnits '^2'],'0-1','Degrees',pixUnits,pixUnits};%Respective units
geoPropConv = {@(x)(x*pixSize^2),@(x)(x), @(x)(x), @(x)(x*pixSize), @(x)(x*pixSize)};%And functions for converting to these units
nGeoProp = numel(geoPropFields);



for j = 1:nGeoProp
    geoProps.(geoPropFields{j}) = nan(nAdhesions,nFrames);
end


for iFrame = 1:nFrames
    
    currMask = segProc.loadChannelOutput(p.MaskChannelIndex,iFrame);
                
    %Get geometric properties of each adhesioin
    currProp = regionprops(currMask,geoPropFields);                           
    
    %Store them in common array.
    currInd = unique(currMask(currMask(:)>0));%Indices present in current frame        
    
    for j = 1:nGeoProp
        
        geoProps.(geoPropFields{j})(currInd,iFrame) = arrayfun(geoPropConv{j},[currProp(currInd).(geoPropFields{j})]);
        
    end    
    
end

varToSave = [varToSave {'geoProps'}];

% ------- Per-Adhesion Figures ------- %

for j = 1:nGeoProp
    
    gFig = fsFigure(.5);    
    imHan = imagesc(tData,1:nAdhesions,geoProps.(geoPropFields{j}));
    gAx = get(gFig,'CurrentAxes');
    saturateImageColormap(gAx,p.SatPercent)
    colormap(p.ColorMap(nCols))
    xlabel(['Time, ' tUnits])
    ylabel('Adhesion #')
    axis xy
    set(imHan,'AlphaData',~isnan(geoProps.(geoPropFields{j})))
    colorbar('location','EastOutside')
    title(['Per-Adhesion ' geoPropFields{j} ' over time, ' geoPropUnits{j}])
    
    if p.SaveFigures; mfFigureExport(gFig,[p.OutputDirectory filesep 'per-adhesions ' geoPropFields{j} ' over time']);end
    
    
end



%% ------------ Temporal Properties ----------- %%

disp('Calculating adhesion temporal properties....')

%Use an arbitrary field to determine during which frames adhesions exist
existMat = ~isnan(geoProps.(geoPropFields{1}));

tPropFields = {'StartTime','EndTime','Lifetime'};
iLife = 3;%Lifetime is special case
tPropFuns = {@(x)(find(x,1,'first')) ,@(x)(find(x,1,'last')), @(x)(nnz(x)*tInterval)};
tPropUnits = {'Frame #','Frame #',tUnits};
ntProp =numel(tPropFields);


for j = 1:nAdhesions
    
    for k = 1:ntProp
        if nnz(existMat(j,:)) > 0
            tProps.(tPropFields{k})(j) = tPropFuns{k}(existMat(j,:));
        end
    end
    
end

varToSave = [varToSave {'tProps'}];

%Make generic figures for all adhesions
for k = 1:ntProp
    
   tFig = figure;
   hist(tProps.(tPropFields{k}),20)
   xlabel([tPropFields{k} ', ' tPropUnits{k}])
   ylabel('# of Adhesions')       
   title('Temporal properties, ALL detected adhesions')
   if p.SaveFigures; mfFigureExport(tFig,[p.OutputDirectory filesep 'all adhesions ' tPropFields{k} ' histogram ']);end
        
end

% ----- Lifetime Figure ----- %

%And one specific figure for the adhesions which started and ended in the
%movie
isComplete = tProps.StartTime > 1 & tProps.EndTime < nFrames;
nComplete = nnz(isComplete);
cFig = figure;
hist(tProps.(tPropFields{iLife})(isComplete),20)
xlabel([tPropFields{iLife} ', ' tPropUnits{iLife}])
ylabel('# of Adhesions') 
title({['Lifetimes, n=' num2str(nComplete) ' complete adhesions'],'(Start and end within movie)'})

if p.SaveFigures, mfFigureExport(cFig,[p.OutputDirectory filesep 'Lifetimes of complete adhesions']);end

% ----- N over Time Figure ----- %

nAdhPerFrame = sum(existMat,1);

nVtFig = figure;
plot(tData,nAdhPerFrame)
xlabel(['Time, ' tUnits])
ylabel('Total # Adhesions')
title({'Adhesion Count over Time',['Total of n=' num2str(nAdhesions) ' unique adhesions']})
if p.SaveFigures, mfFigureExport(nVtFig,[p.OutputDirectory filesep 'Adhesion count vs time']);end


%% ------------- Alignment --------- %%
%Aligns all complete adhesions to start time


% ---- Aligned, complete adhesions ---- %
iComplete = find(isComplete);
for j= 1:nGeoProp    
    alignedGeoProps.(geoPropFields{j}) = nan(nComplete,nFrames);
end
for j= 1:nIntProp
    for l = 1:nChan
        alignedIntProps(l).(intPropFields{j}) = nan(nComplete,nFrames);
    end
end

alignedN = nan(nComplete,nFrames);

for i = 1:nComplete
    
    currTind = tProps.StartTime(iComplete(i)):tProps.EndTime(iComplete(i));
    currI = iComplete(i);
    currAlignT = 1:round(tProps.Lifetime(currI)/tInterval);
    
    for j = 1:nGeoProp                                  
        alignedGeoProps.(geoPropFields{j})(i,currAlignT) = geoProps.(geoPropFields{j})(currI,currTind);                    
    end
    
    for l = 1:nChan
        %TEMP - convert to using function handles!!?
        alignedIntProps(l).MeanIntensity(i,currAlignT) = squeeze(allSamples(l).avg(currI,1,currTind));            
    end
    
    %Get aligned # of pixels for weighting of combined means (same for both
    %channels)
    alignedN(i,currAlignT) = squeeze(allSamples(1).n(currI,1,currTind));
    
end

varToSave = [varToSave {'alignedGeoProps','alignedIntProps','alignedN'}];

% ---- Aligned, per-complete-adhesion figures --- %

for j = 1:nGeoProp
    
    gFig = fsFigure(.5);    
    imHan = imagesc(tData,1:nComplete,alignedGeoProps.(geoPropFields{j}));
    gAx = get(gFig,'CurrentAxes');
    saturateImageColormap(gAx,p.SatPercent)
    colormap(p.ColorMap(nCols))
    xlabel(['Time, ' tUnits])
    ylabel('Complete Adhesion #')
    axis xy
    set(imHan,'AlphaData',~isnan(alignedGeoProps.(geoPropFields{j})))
    colorbar('location','EastOutside')
    title(['Aligned, Per-Complete-Adhesion ' geoPropFields{j} ' over time, ' geoPropUnits{j}])
    
    if p.SaveFigures, mfFigureExport(gFig,[p.OutputDirectory filesep 'Aligned per-complete adhesion ' geoPropFields{j} ' vs time']);end
    
end

for j = 1:nIntProp
    
    for k=  1:nChan
        gFig = fsFigure(.5);    
        imHan = imagesc(tData,1:nComplete,alignedIntProps(k).(intPropFields{j}));
        gAx = get(gFig,'CurrentAxes');
        saturateImageColormap(gAx,p.SatPercent)
        colormap(p.ColorMap(nCols))
        xlabel(['Time, ' tUnits])
        ylabel('Complete Adhesion #')
        axis xy
        set(imHan,'AlphaData',~isnan(alignedIntProps(k).(intPropFields{j})))
        colorbar('location','EastOutside')
        title(['Aligned, Per-Complete-Adhesion ' intPropFields{j} ' in channel ' num2str(p.ChannelIndex(k)) ' over time, ' intPropUnits{j}])
        
        if p.SaveFigures, mfFigureExport(gFig,[p.OutputDirectory filesep 'Aligned per-complete adhesion ' intPropFields{j} ' vs time']);end
    end
end


%% ----- Combined Geometric & Intensity Properties over Time ------- %%



% ---- Make figures, calc Unaligned Stats---- %

combTPropFields = {'Mean','STD','Total'};
combTPropFun = {@(x,y)(nansum( x .* y,1) ./ nansum(y,1)), @(x,y)(sqrt(nanvarXT(x,y,1))), @(x,y)(nansum(x .* y,1))};%x - values, y - weights (# pixels)
nCombTProp = numel(combTPropFields);

for j = 1:nGeoProp
    
    %Get weights for current field
    if ~strcmp(geoPropFields{j},'Area')            
        weightsAll = squeeze(allSamples(1).n(:,1,:));
        weightsAlign = alignedN;            
        wtStr = 'area-weighted, ';
    else
        %We don't weight the area by the area.....
        weightsAll = ones(nAdhesions,nFrames);
        weightsAlign = ones(nComplete,nFrames);
        wtStr = '';
    end
       
    for k = 1:nCombTProp
        combTPropNames{j,k} = [combTPropFields{k} geoPropFields{j}];
        
        combTProps.(combTPropNames{j,k}) = combTPropFun{k}(geoProps.(geoPropFields{j}),weightsAll);                        
        combTPropsAligned.(combTPropNames{j,k}) = combTPropFun{k}(alignedGeoProps.(geoPropFields{j}),weightsAlign);
    end
    
    ctFig = fsFigure(.6);
    hold on    
    
    subplot(2,1,1)
    hold on   
    title({['Combined, unaligned, ' wtStr 'mean of ' geoPropFields{j} ', All adhesions, Total of n=' num2str(nAdhesions) ' unique adhesions']})            
    plot(tData,combTProps.(combTPropNames{j,1}))
    plot(tData,combTProps.(combTPropNames{j,1}) + combTProps.(combTPropNames{j,2}),'--')
    plot(tData,combTProps.(combTPropNames{j,1}) - combTProps.(combTPropNames{j,2}),'--')
    legend('Mean','+/-STD')
    xlabel(['Time, ' tUnits])
    ylabel([combTPropNames{j,1} ', ' geoPropUnits{j}])        

    subplot(2,1,2)
    hold on   
    title({['Combined, unaligned, ' wtStr ' sum of ' geoPropFields{j} ', All adhesions, Total of n=' num2str(nAdhesions) ' unique adhesions']})            
    plot(tData,combTProps.(combTPropNames{j,3}))        
    xlabel(['Time, ' tUnits])
    ylabel([combTPropNames{j,3} ', ' geoPropUnits{j}])        
    
    if p.SaveFigures, mfFigureExport(ctFig,[p.OutputDirectory filesep 'Combined '  geoPropFields{j} ' of all adhesions vs time']);end
    
    ctaFig = fsFigure(.6);
    hold on       
    
    subplot(2,1,1)
    hold on   
    title({['Combined, aligned, ' wtStr 'mean of ' geoPropFields{j} ', complete adhesions, Total of n=' num2str(nComplete) ' complete adhesions']})            
    plot(tData,combTPropsAligned.(combTPropNames{j,1}))
    plot(tData,combTPropsAligned.(combTPropNames{j,1}) + combTPropsAligned.(combTPropNames{j,2}),'--')
    plot(tData,combTPropsAligned.(combTPropNames{j,1}) - combTPropsAligned.(combTPropNames{j,2}),'--')
    legend('Mean','+/-STD')
    xlabel(['Time, ' tUnits])
    ylabel([combTPropNames{j,1} ', ' geoPropUnits{j}])        

    subplot(2,1,2)
    hold on   
    title({['Combined, aligned, ' wtStr 'sum of ' geoPropFields{j} ', complete adhesions, Total of n=' num2str(nComplete) ' complete adhesions']})            
    plot(tData,combTPropsAligned.(combTPropNames{j,3}))        
    xlabel(['Time, ' tUnits])
    ylabel([combTPropNames{j,3} ', ' geoPropUnits{j}])        
    
    if p.SaveFigures, mfFigureExport(ctaFig,[p.OutputDirectory filesep 'Combined aligned '  geoPropFields{j} ' of complete adhesions vs time']);end
        
    
end

varToSave = [varToSave {'combTProps'}];

for iChan = 1:nChan
    for j = 1:nIntProp

        weightsAlign = alignedN;            
        wtStr = 'area-weighted, ';    


        for k = 1:nCombTProp

            combIntPropNames{j,k} = [combTPropFields{k} intPropFields{j}];
            combIntPropsAligned(iChan).(combIntPropNames{j,k}) = combTPropFun{k}(alignedIntProps(iChan).(intPropFields{j}),weightsAlign);        

        end 

        ctaiFig = fsFigure(.6);
        hold on   

        subplot(2,1,1)
        hold on   
        title({['Channel ' num2str(p.ChannelIndex(iChan)) ' Combined, aligned, ' wtStr 'mean of ' intPropFields{j} ', complete adhesions, Total of n=' num2str(nComplete) ' complete adhesions']})            
        plot(tData,combIntPropsAligned(iChan).(combIntPropNames{j,1}))
        plot(tData,combIntPropsAligned(iChan).(combIntPropNames{j,1}) + combIntPropsAligned(iChan).(combIntPropNames{j,2}),'--')
        plot(tData,combIntPropsAligned(iChan).(combIntPropNames{j,1}) - combIntPropsAligned(iChan).(combIntPropNames{j,2}),'--')
        legend('Mean','+/-STD')
        xlabel(['Time, ' tUnits])
        ylabel([combIntPropNames{j,1} ', ' intPropUnits{j}])        

        subplot(2,1,2)
        hold on   
        title({['Channel ' num2str(p.ChannelIndex(iChan)) ' Combined, aligned, ' wtStr 'sum of ' intPropFields{j} ', complete adhesions, Total of n=' num2str(nComplete) ' complete adhesions']})            
        plot(tData,combIntPropsAligned(iChan).(combIntPropNames{j,3}))        
        xlabel(['Time, ' tUnits])
        ylabel([combIntPropNames{j,3} ', ' intPropUnits{j}])        

        if p.SaveFigures, mfFigureExport(ctaiFig,[p.OutputDirectory filesep 'Channel ' num2str(p.ChannelIndex(iChan)) ' Combined aligned '  intPropFields{j} ' of complete adhesions vs time']);end

    end
end

varToSave = [varToSave {'combIntPropsAligned'}];

%% ----------------------- Export to Spreadsheet --------------------- %%

geoToExport = {'Area'};
intToExport = {'MeanIntensity'};
fExt = '.csv';

for j = 1:numel(geoToExport)

    %Export raw
    iExp = find(strcmp(geoToExport{j},geoPropFields));
    expFileName = [p.OutputDirectory filesep 'per adhesion ' geoPropFields{iExp} ' ' geoPropUnits{iExp} fExt];    
    outArray = geoProps.(geoPropFields{iExp});    
    csvwrite(expFileName,outArray);
        
    %And aligned
    expFileName = [p.OutputDirectory filesep 'per aligned complete adhesion ' geoPropFields{iExp} ' ' geoPropUnits{iExp} fExt];    
    outArray = [iComplete(:)  alignedGeoProps.(geoPropFields{iExp})];    
    csvwrite(expFileName,outArray);
        
end


for j = 1:numel(intToExport)
    
    for k = 1:nChan

        iExp = find(strcmp(intToExport{j},intPropFields));
        expFileName = [p.OutputDirectory filesep 'per aligned complete adhesion ' intPropFields{iExp} ' channel ' num2str(p.ChannelIndex(k))  fExt];

        outArray = [iComplete(:) alignedIntProps(k).(intPropFields{iExp})];

        csvwrite(expFileName,outArray);
    end
end

%And write the raw intensity samples as well
for j = 1:nChan
            
    expFileName = [p.OutputDirectory filesep 'per adhesion ' intPropFields{1} ' channel ' num2str(p.ChannelIndex(j)) fExt];
    outArray = squeeze(allSamples(j).avg);
    
    csvwrite(expFileName,outArray);
    
end

%% ------ Export to .mat ------ %%


save([p.OutputDirectory filesep outFileName '.mat'],varToSave{:})

disp('Finished adhesion analysis!')
disp(['Results written to : ' p.OutputDirectory])

