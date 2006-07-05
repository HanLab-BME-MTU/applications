function [idlist,dataPropertiesOut,success] = LG_recalcIdlist(idlist,dataProperties,recalcOptions)
%LG_recalcIdlist will call linker to redo idlist


% init success, dataPropertiesOut
success = 0;
dataPropertiesOut = [];


% define order of recalcOptions (fieldnames in dataProperties)
recalcNames = cell(8,1);
recalcNames{1} = 'recalc';
recalcNames{3} = 'MAXSPOTS';
recalcNames{4} = 'linker_relAmpWeight';
recalcNames{5} = 'linker_useCOM';
recalcNames{6} = 'linker_relativeMaxDistance';
recalcNames{7} = 'linker_absoluteMaxDistance';
recalcNames{2} = 'fuseTags';
recalcNames{8} = 'linker_fuseRatio';


% initialize dataPropertiesHasChanged
dataPropertiesHasChanged = 0;

if nargin < 3 || isempty(recalcOptions)
    
    % add (new) default fields to dataProperties
    dataPropertiesNew = defaultDataProperties(dataProperties);
    if ~isequalwithequalnans(dataPropertiesNew,dataProperties)
        dataProperties = dataPropertiesNew;
        dataPropertiesHasChanged = 1;
    end
    clear dataPropertiesNew

    % ask for options - make a nice gui later
    prompt={['Select the extent of the recalc',char(10),...
        '0: Full recalc',char(10),...
        '1: Keep maxSpot-frame-links',char(10),...
        '2: Only estimate positions and amplitudes',char(10)],...
        'fuse tags',...
        'MaxNumSpots',...
        sprintf('Weight of amplitude relative to distance (currently 1/%1.3f)',...
        dataProperties.linker_relAmpWeight),...
        'Use COM to correct positions',...
        'relative distance cutoff',...
        'absolute distance cutoff (\mum)',...
        'ratio of Rayleigh-limit by which the tags have to be separated to still allow fusion'};
    name='Set recalc options';
    numlines=1;
    defaultAnswer={'0','0',num2str(dataProperties.MAXSPOTS),...
        sprintf('%1.3f',dataProperties.linker_relAmpWeight),...
        sprintf('%i',dataProperties.linker_useCOM),...
        sprintf('%i',dataProperties.linker_relativeMaxDistance),...
        sprintf('%i',dataProperties.linker_absoluteMaxDistance),...
        sprintf('%1.2f',dataProperties.linker_fuseRatio)};

    options.Interpreter='tex';

    answer=inputdlg(prompt,name,numlines,defaultAnswer,options);

    % quit if cancel
    if isempty(answer)
        return
    end

    % loop through answers. check for empty
    emptyAnswer = 0;
    for i = length(defaultAnswer):-1:1
        answerI = str2double(answer{i});
        if isempty(answerI) || isnan(answerI)
            emptyAnswer = 1;
        else
            if i>1
            % check whether there has been a change at all
            if isequalwithequalnans(defaultAnswer(i),answer(i))
                recalcOptions{i} = NaN;
            else
                recalcOptions{i} = answerI;
            end
            else
                % don't check for change in 'recalc'!
                recalcOptions{i} = answerI;
            end
        end
    end

    if emptyAnswer
        h = errordlg('Empty input field or non-numeric input!','Error processing recalcOptions');
        uiwait(h);
        return
    end
    
    % with the dialogue, we will always have the full complement of options
    nOptions = 8;
    
else
    
    % count number of options
    nOptions = length(recalcOptions);
    if nOptions > 8
        h = errordlg(sprintf('Too many recalc options (%i instead of 8)',nOptions),...
            'Error processing recalcOptions');
        uiwait(h);
        return
    end
end

% CHANGE - DATAPROPERTIES SHOULD ALWAYS REFLECT LAST GOOD SET OF OPTIONS!

% loop through recalcNames. write recalcOptions into dataProperties, and
% keep recalcOptions{1}, which lists what exactly to recalc, and keep
% recalcOptions{2}, which tells whether to fuse or not

for r=nOptions:-1:3
    % only add option if it's not NaN
    if ~isnan(recalcOptions{r})
        dataProperties.(recalcNames{r}) = recalcOptions{r};
        % remove option from the list
        recalcOptions(r)=[];
        % remember that dataProperties have changed
        dataPropertiesHasChanged = 1;
    else
        % just remove entry from list
        recalcOptions(r) = [];
    end
end

% add 10 to recalcOptions{1} in case we fuse tags
if length(recalcOptions) > 1
    recalcOptions{1} = recalcOptions{1} + 10 * (recalcOptions{2}>0);
end

if ~iscell(recalcOptions)
    recalcOptions = {recalcOptions};
end

% now that all is hopefully good, recalcIdlist
idlist(1).stats.recalc = recalcOptions(1);
idlist = linker(idlist, dataProperties);

if ~dataPropertiesHasChanged
    dataPropertiesOut = [];
else
    dataPropertiesOut = dataProperties;
end

success = 1;