function [idlist,dataProperties,success] = LG_recalcIdlist(idlist,dataProperties,recalcOptions)
%LG_recalcIdlist will call linker to redo idlist

% define order of recalcOptions (fieldnames in dataProperties)
recalcNames = cell(6,1);
recalcNames{1} = 'recalc';
recalcNames{2} = 'MAXSPOTS';
recalcNames{3} = 'linker_relAmpWeight';
recalcNames{4} = 'linker_useCOM';
recalcNames{5} = 'linker_relativeMaxDistance';
recalcNames{6} = 'linker_absoluteMaxDistance';


% initialize dataPropertiesHasChanged
dataPropertiesHasChanged = 0;

if nargin < 3 || isempty(recalcOptions)
    
    % add (new) default fields to dataProperties
    dataPropertiesNew = defaultDataProperties(dataProperties);
    if ~isequalwithequalnans(dataPropertiesNew,dataProperties)
        dataProperties = dataPropertiesNew;
        dataPropertiesHasChanged = 1;
    end
    dataPropertiesNew = [];

    % ask for options - make a nice gui later
    prompt={['Select the extent of the recalc',char(10),...
        '0: Full recalc',char(10),...
        '1: Keep maxSpot-frame-links',char(10),...
        '2: Only estimate positions and amplitudes',char(10)],...
        'MaxNumSpots',...
        sprintf('Weight of amplitude relative to distance (currently 1/%1.3f)',...
        dataProperties.linker_relAmpWeight),...
        'Use COM to correct positions',...
        'relative distance cutoff',...
        'absolute distance cutoff (um)'};
    name='Set recalc options';
    numlines=1;
    defaultAnswer={'0',num2str(dataProperties.MAXSPOTS),...
        sprintf('%1.3f',dataProperties.linker_relAmpWeight),...
        sprintf('%i',dataProperties.linker_useCOM),...
        sprintf('%i',dataProperties.linker_relativeMaxDistance),...
        sprintf('%i',dataProperties.linker_absoluteMaxDistance)};

    options.Interpreter='tex';

    answer=inputdlg(prompt,name,numlines,defaultAnswer,options);

    % quit if cancel
    if isempty(answer)
        success = 0;
        dataProperties = [];
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
        success = 0;
        dataProperties = [];
        return
    end
    
    % with the dialogue, we will always have the full complement of options
    nOptions = 6;
    
else
    
    % count number of options
    nOptions = length(recalcOptions);
    if nOptions > 6
        h = errordlg(sprintf('Too many recalc options (%i instead of 6)',nOptions),...
            'Error processing recalcOptions');
        uiwait(h);
        success = 0;
        dataProperties = [];
        return
    end
end

% CHANGE - DATAPROPERTIES SHOULD ALWAYS REFLECT LAST GOOD SET OF OPTIONS!

% loop through recalcNames. write recalcOptions into dataProperties, and
% keep only recalcOptions{1}, which lists what exactly to recalc

for r=nOptions:-1:2
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

% now that all is hopefully good, recalcIdlist
idlist(1).stats.recalc = recalcOptions(1);
idlist = linker(idlist, dataProperties);

if ~dataPropertiesHasChanged
    dataProperties = [];
end

success = 1;