function [ output_args ] = GCAAnalysisExtract_VeilStemThicknessMovie(movieData,varargin)
% GCAAnalysisExtract_VeilStemThickness
%%  small wrapper to extract the veil thickness


%% PARAMS
%
% %% PARAMS: %%
%
%    'distCutOffForThickness' (PARAM) : Positive scalar
%        The distance along the longest path in the veil/stem skeleton
%        from the tip of the leading protrusion to the this user defined
%        point.
%        Default : 20  (currently in um)
%
%% CheckInput
ip = inputParser;

ip.CaseSensitive = false;
ip.addOptional('veilStem',[]); % have option to load
ip.addOptional('neuriteLength',[]);
ip.addParameter('distCutOffForThickness',20); % in um
ip.addParameter('visual',true);

ip.parse(varargin{:});
p = ip.Results;
pToSave = rmfield(p,{'veilStem','neuriteLength'});
%%
veilStem = ip.Results.veilStem{1};
neuriteLength = ip.Results.neuriteLength{1};
%% FIX  load the neurite length values and veil stem if not input
% if isempty(veilStem)....
%%





%% Check if Step IV run





%%
nFrames = movieData.nFrames_;
measC = cell(nFrames,1);
outDir = [movieData.outputDirectory_ filesep ...
    'MEASUREMENT_EXTRACTION' filesep 'Descriptor' filesep 'Veil' filesep ...
    'VeilStemThickness'];

if ~isempty(outDir)
    mkdir(outDir)
end

idx  = neuriteLength > ip.Results.distCutOffForThickness;
if (sum(idx) == length(neuriteLength));
     
else
    newVal = max(neuriteLength); 
   pToSave.distCutOffForThickness  = max(neuriteLength);
    errorMessage =[ 'Max Neurite Length Sampled Less than Value for distCutOffForThickness: Resetting Value to ' num2str(newVal)] ;
    display(errorMessage);
    save([outDir filesep 'ErrorReport.mat'],'errorMessage');
end

if ip.Results.visual == true;
    measDirVis = [outDir filesep 'MeasurementVisual']; 
    if ~isdir(measDirVis)
        mkdir(measDirVis); 
    end    
end


for iFrame = 1:(nFrames-1)
    longPathLinIndC=  veilStem(iFrame).neuriteLongPathIndices;
    veilStemMaskC = veilStem(iFrame).finalMask;
    % get the thickness metric
    [thicknessValues,TSFig]=   GCAAnalysisExtract_VeilStemThickness(longPathLinIndC,veilStemMaskC,pToSave);
    % put into format
    measC{iFrame,1} = thicknessValues;
    if ip.Results.visual == true
        
        saveas(TSFig.h,[measDirVis num2str(iFrame,'%03d') '.png' ]);
         
    end
    clear TSFig
end

save([outDir filesep 'meas_VeilStemThickness.mat'],'measC');
% save the input
save([outDir filesep 'p.mat'],'pToSave');
end





