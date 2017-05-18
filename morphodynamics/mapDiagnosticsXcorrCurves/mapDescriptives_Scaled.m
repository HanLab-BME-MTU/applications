function mapDescriptives_Scaled(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, varargin)
% mapDescriptives_Scaled Draw descriptive plots of an scaled activity map of
% the specified channel in movieData.
%
%
%
% Jungsik Noh, 2017/01/09


ip = inputParser;
%ip.addParameter('adf', 0);
ip.addParameter('figFlag', 'off');
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
%ip.addParameter('parpoolNum', 4);
%ip.addParameter('rseed', 'shuffle');
%ip.addParameter('numPerm', 1000);
ip.addParameter('omittedWindows', []);

ip.parse(varargin{:});
p = ip.Results;

figFlag = p.figFlag;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isdir(figuresDir); mkdir(figuresDir); end


%%  getting Maps from channels

disp(chanName)
disp(chanTitle)


[fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, ~] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
            'WithN', p.WithN, 'omittedWindows', p.omittedWindows); 
disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])


%%  scaling the activity map

actmap_outlSc = cell(1, maxLayer);
for indL = 1:maxLayer
    
    actmap_outlCell = num2cell(actmap_outl{indL}, 2);
    actmap_outlZ = cellfun(@(x) nanZscore(x), actmap_outlCell, 'UniformOutput', false);
    actmap_outlSc{indL} = cell2mat(actmap_outlZ);
    
    % CentMap ...
    %actmap_outlSc{indL} = detrend(actmap_outl{indL}', 'constant')';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plots 

%%  smoothActivityMap prot/act maps

smParam = 0.8;

fchan = cell(1, maxLayer);
for indL = 1:maxLayer
    
    inputmap = actmap_outlSc{indL};

    filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
    fchan{indL} = figure('Visible', figFlag);
    figtmp = imagesc(filteredmap);
    title([chanTitle, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)

    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;

end


%% 
for indL = 1:maxLayer
    saveas(fchan{indL}, fullfile(figuresDir, [fname0, 'ScMap_', num2str(indL), 'L.png']), 'png')
end




end
