function [ h ] = GCAAnalysisInitialReport(projPath,toPlotGroup,iProj)
%display(strrep(toPlotGroup.info.projList{1}{iProj,2},'_',''));

%params{1} = 'expression';
params{1} = 'filopodia';
params{2} = 'density';





% options: 'veil', 'density', 'orientation',expression

load([projPath filesep 'ANALYSIS' filesep 'movieData.mat']);
% fonts for distribution plots within group
axisFont = 14;
labelFont = 20;
titleFont = 30;
[~,Exp] = upDirectory(projPath,2,1);

%load([projPath filesep 'ANALYSIS' filesep 'notes.mat']);
% load all the notes
notesAll = cellfun(@(x) load([x filesep 'ANALYSIS' filesep 'notes.mat']),toPlotGroup.info.projList{1}(:,1));
% filter based on the current include tag
idxInclude =   arrayfun(@(x) ~isempty(regexpi(notesAll(x).notes.Include_Currently,'Yes')),1:length(notesAll));
toPlotGroup.info.projList{1} = toPlotGroup.info.projList{1}(idxInclude,:);
toPlotGroup.info.grouping = ones(sum(idxInclude),1);
projListC = toPlotGroup.info.projList{1}(:,1);

idxC = find(cellfun(@(x) strcmpi(x,projPath),projListC));
notes = notesAll(iProj).notes;
%% Cropping

fsFigure(1);

% load cropped overlay
subplot(1,2,1);
%img1 = imread([projPath filesep 'CropOverlay' filesep '001.tif']);
%imshow(img1);
% openfig([projListC{iProj} filesep 'CropOverlay' filesep '001.fig'],'new','visible');

plotCroppedRegion(projPath);


title({'Cropping Frame Overlay:' ; 'Frame 1'});
pixels = round(10/0.216); % make a variable both scale bar in micron and the pixel size: scale bar set to 10 um here
plotScaleBar(pixels,pixels/20,'Color',[0 0 0],'Location', 'SouthWest');

subplot(1,2,2);
img =  double(imread(([projPath filesep 'Channels' filesep 'C1_mCherry' filesep ...
    'C1_001.tif'])));
imshow(-img,[]);
pixels = round(10/0.216); % make a variable both scale bar in micron and the pixel size: scale bar set to 10 um here
plotScaleBar(pixels,pixels/20,'Color',[0 0 0],'Location', 'SouthWest');
title({'Raw Cropped Image:' ;'Frame 1'});
set(gcf,'Color',[1,1,1]); % make sure it is set to white.
snapnow
close gcf a
%% Segmentation
if ~isempty(regexpi(notes.Include_Currently,'Yes'));
    toPlotGroup = GCACollectGroupData(toPlotGroup); 
    toPlotGroup = GCAReformatVeilStats(toPlotGroup,0);
    fsFigure(1);
    subplot(1,2,1);
    %img3 = imread([projPath filesep 'ANALYSIS' filesep 'CompareOtsuVeilStem' ...
    %    filesep  'Compare001.png']);
    %openfig([projPath filesep 'ANALYSIS' filesep 'CompareOtsuVeilStem' ...
    % 'Compare001.fig']);
    wrapCompareGlobalOtsuVeilStem(MD,1);
    title({'Comparison to Global Thresholding:' ; 'Frame 1'});
    
    subplot(1,2,2);
    
    img4 = imread([projPath filesep 'ANALYSIS' filesep 'Visualization_Overlays'...
        filesep 'Overlaysmono' filesep '001.png']);
    % openfig([projPath filesep 'ANALYSIS' filesep 'growth_cone_mask'...
    %   'Overlaysmono' filesep '001.fig']);
    imshow(img4);
    title({'Segmentation Overlay:' ; 'Frame 1'});
    set(gcf,'Color',[1,1,1]); % make sure it is set to white.
    
    snapnow
    close gcf
end
%% NOTES MB
%%   * Morphology
if ~isempty(notes.Morphology);
    if length(notes.Morphology)>60;
        num = ceil(length(notes.Morphology)/60);
        start = 1;
        % find spaces
        spaces = regexp(notes.Morphology,' ');
        for i = 1:num
            
            
            %start = find((spaces>start & spaces<50*num),1,'first');
            stop = find(spaces<60*i,1,'last');
            if stop == length(spaces);
                finish = length(notes.Morphology);
            else
                finish = spaces(stop);
            end
            display(notes.Morphology(start:finish));
            start = finish+1;
            
        end
    else
        display(notes.Morphology)
    end
end
%%   * Crop
if ~isempty(notes.Crop);
    if  length(notes.Crop)>60
        
        
        num = ceil(length(notes.Crop)/60);
        start = 1;
        % find spaces
        spaces = regexp(notes.Crop,' ');
        for i = 1:num
            
            
            %start = find((spaces>start & spaces<50*num),1,'first');
            stop = find(spaces<60*i,1,'last');
            if stop == length(spaces);
                finish = length(notes.Crop);
            else
                finish = spaces(stop);
            end
            display(notes.Crop(start:finish));
            start = finish+1;
            
        end
    else
        display(notes.Crop)
    end
end
%%   * Veil/Stem
if ~isempty(notes.VeilStem);
    if length(notes.VeilStem) > 60
        num = ceil(length(notes.VeilStem)/60);
        start = 1;
        % find spaces
        spaces = regexp(notes.VeilStem,' ');
        for i = 1:num
            
            
            %start = find((spaces>start & spaces<50*num),1,'first');
            stop = find(spaces<60*i,1,'last');
            if stop == length(spaces);
                finish = length(notes.VeilStem);
            else
                finish = spaces(stop);
            end
            display(notes.VeilStem(start:finish));
            start = finish+1;
            
        end
    else
        display(notes.VeilStem);
    end
end
%%   * Filopodia
if ~isempty(notes.Filopodia);
    if length(notes.Filopodia)>60;
        num = ceil(length(notes.Filopodia)/60);
        start = 1;
        % find spaces
        spaces = regexp(notes.Filopodia,' ');
        for i = 1:num
            
            
            %start = find((spaces>start & spaces<60*num),1,'first');
            stop = find(spaces<60*i,1,'last');
            if stop == length(spaces);
                finish = length(notes.Filopodia);
            else
                finish = spaces(stop);
            end
            display(notes.Filopodia(start:finish));
            start = finish+1;
            
        end
    else
        display(notes.Filopodia);
    end
end
%%   * Include Currently
if ~isempty(notes.Include_Currently);
    if length(notes.Include_Currently)>60;
        num = ceil(length(notes.Include_Currently)/60);
        start = 1;
        % find spaces
        spaces = regexp(notes.Include_Currently,' ');
        for i = 1:num
            
            
            %start = find((spaces>start & spaces<60*num),1,'first');
            stop = find(spaces<60*i,1,'last');
            if stop == length(spaces);
                finish = length(notes.Include_Currently);
            else
                finish = spaces(stop);
            end
            display(notes.Include_Currently(start:finish));
            start = finish+1;
            
        end
    else
        display(notes.Include_Currently);
    end
end
%%   * Specific Questions
if ~isempty(notes.Specific_Questions);
    if length(notes.Specific_Questions)>60;
        
        num = ceil(length(notes.Specific_Questions)/60);
        start = 1;
        % find spaces
        spaces = regexp(notes.Specific_Questions,' ');
        for i = 1:num
            
            
            %start = find((spaces>start & spaces<60*num),1,'first');
            stop = find(spaces<60*i,1,'last');
            if stop == length(spaces);
                finish = length(notes.Specific_Questions);
            else
                finish = spaces(stop);
            end
            display(notes.Specific_Questions(start:finish));
            start = finish+1;
            
        end
    else
        display(notes.Specific_Questions);
    end
end


%% NOTES PERTZ - Please feel free to add any comments you may have using the pdf browser. Thanks!
% * Morphology?
% * Crop?
% * Veil/Stem?
% * Filopodia?
% * Exclude Final?
if ~isempty(regexpi(notes.Include_Currently,'Yes'));
    
    
    %% Visual Outgrowth : Veil : Filopodia In Time
    
    fsFigure(1);
    load([projPath filesep 'ANALYSIS' filesep 'PARAMETER_EXTRACTION' ...
        filesep 'GlobalFunctional\neurite_outgrowth_measurements' filesep 'neuriteLengthOutput.mat'])
    subplot(3,1,1);
    
    plotNeuriteOutgrowthInTime(neuriteLength,'k',1,5,0,[],'subplot',0);
    title('Neurite Outgrowth','FontName','Arial','FontSize',12);
    
    subplot(3,1,2);
    if exist([projPath filesep 'ANALYSIS' filesep 'protrusion_samples' filesep ...
            'protrusion_samples.mat'],'file')~=0;
        
        load([projPath filesep 'ANALYSIS' filesep 'protrusion_samples' filesep ...
            'protrusion_samples.mat']);
        protrusionMapMine(protSamples,[],0);
        title('Local Veil Dynamics','FontName','Arial','FontSize',12);
    end
    subplot(3,1,3);
    % img = imread([projPath filesep 'ANALYSIS\PARAMETER_EXTRACTION\Filopodia\Length\External' filesep '001.png']);
    % imshow(img);
    load([projPath filesep 'ANALYSIS\PARAMETER_EXTRACTION\Descriptor\Filopodia\ConnectToVeil_LengthInt\filoLength' filesep ...
        'param_filoLengthToVeil.mat']);
    distribInTime = reformatDataCell(paramC); 
    makeMovieBoxplot(distribInTime,'Filopodia Length (um)');
    title('Filopodia Distribution In Time','FontName','Arial','FontSize',12)
    snapnow
    close gcf
%     %% Quantification With Respect to Group : Check Fluorescent Protein Expression Level
%     fsFigure(1);
%     projListAll = vertcat(toPlotGroup.info.projList{:});
%     
%     % Filopodia Orientation Relative to Local Veil (Degrees)
%     toPlotGroup = reformatExpressionLevelForPlotting(toPlotGroup,0);
%     names = projListAll(:,2);
%     names = strrep(names,'_',' ');
%     dataMat = toPlotGroup.expressLevel{1};
%     h12 = boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
%     %h1=  boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
%     %set(h1(:),'Linewidth',1);
%     values = dataMat(:);
%     values = values(~isnan(values));
%     minValue = min(values);
%     maxValue = max(values);
%     
%     
%     set(h12(:),'Linewidth',1);
%     idx = 1:(size(dataMat,2));
%     idx = idx==idxC;
%     dataMat(1:size(dataMat,1),~idx)=NaN;
%     hold on
%     h22 = boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
%     
%     axis([0.5,size(dataMat,2)+0.5, minValue-5, maxValue+5]);
%     set(h22(:),'Linewidth',5);
%     set(gca,'FontSize',axisFont);
%     title(['Group: ' toPlotGroup.info.names{1} ' ' Exp ],'FontName','Arial','FontSize',titleFont);
%     ylabel({'Mean Background Subtracted Intensity' ; 'of Growth Cone (AU)'}','FontName','Arial','FontSize',labelFont);
%     
%     
%     snapnow
%     
%     close gcf
    
    %% Quantification With Respect to Group : Outgrowth
    
    fsFigure(1);
    subplot(1,2,1);
    
    [neuriteLengths,deltas] = gcaCollectMultNeuriteLengthsAndDeltas(projListC);
    
    hold on
    cellfun(@(x)  plotNeuriteOutgrowthInTime(x,'k',1,5,0,[],[],1,1),neuriteLengths);
    
    % plot thick
    plotNeuriteOutgrowthInTime(neuriteLengths{idxC},'k',1,5,0,[],[],1,5);
    
    
    subplot(1,2,2);
    scatter(ones(length(deltas),1),deltas,100 ,'k','filled','MarkerEdge','w');
    hold on
    scatter(1,deltas(idxC),300,'k','filled','MarkerEdge','w');
    axis([0.5,1.5,-10,25]);
    ylabel({'Delta Neurite Outgrowth (um)' ; 'in 10 mins' },'FontSize',20,'FontName','Arial');
    
    set(gca,'XTick',1);
    set(gca,'XTickLabel',[ toPlotGroup.info.names{1}  ]);
    set(gca,'FontName','Arial','FontSize',18);
    snapnow
    close gcf
    
    %% Quantification With Respect to Group- Filopodia (Traditional Definition: 'External') : Length
   % fsFigure(1);
    %toPlotGroupFilo = reformatLengthForPlotting(toPlotGroup,0); % collect the data per movie ...
    % from each project in the toPlot group and reformat it into the dataMat
    % framework useful for plotting per cell and pooled boxplots.
    fsFigure(1);
    
    projListAll = vertcat(toPlotGroup.info.projList{:}) ;
    names = projListAll(:,2);
    names = strrep(names,'_',' ');
    dataMat = toPlotGroup.filoLengthToVeil{1};
    h1=  boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
    set(h1(:),'Linewidth',1);
    idx = 1:(size(dataMat,2));
    idx = idx==idxC;
    dataMat(1:size(dataMat,1),~idx)=NaN;
    hold on
    h2 = boxplot(dataMat,'color','k','colorGroup',toPlotGroup.info.grouping,'notch','on','outlierSize',1,'Labels',names,'labelorientation','inline');
    set(h2(:),'Linewidth',5);
    set(gca,'FontName','Arial','FontSize',axisFont);
    ylabel('Filopodia Length (um)','FontName','Arial','FontSize',labelFont);
    title(['Group: ' toPlotGroup.info.names{1} ' ' Exp],'FontName','Arial','FontSize',titleFont);
    axis([0.5,length(projListC)+0.5,0,15]);
    snapnow
    close gcf
    %% Quantification With Respect to Group: Veil Protrusion Velocity
    % This Quantification is not yet re-run for this set:
    %
    %
    
end



