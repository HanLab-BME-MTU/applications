function [ toPlot,removeNoCortical,noViableData ] = mitoticGroupAnalysisCollectData(groupList,varargin)
% mitoticGroupAnalysisCollectData: a group analysis function that collects 
% data in a groupList for desired measurement, and stores plotting/grouping information. 
% (ie toPlot structure)
% 
% Currently supports measurements: 
%   percent lateral transitions
%   dwell calcs (shrinkage or terminal event only)
%
% This framework generically feeds into 
% mitoticGroupAnalysisGroupPlots.m so one can make perCell and/or pooled 
% distributions plots/ perform stats for any measurement collected in this 
% format.  
%
% INPUT: 
% groupList: (REQUIRED) : 
% 
% measToCollect: (PARAM) : character ('percentLateral' or 'dwellTime')
%                          Default : 'percentLateral'
%                          
% 
% PARAMETERS FOR CLASS CUTOFFS %
% minDispVect: (PARAM) : numeric
%                        Default : 3 pixels
%                        min MT Track Length
%                        (currently in pixels)
%                        see Supp Figure 1 A of Kwon et al. Dev Cell 2015
%
% orientCutOff: (PARAM) : numeric 
%                         Default : 60 (in degrees) 
%                         Local Direction MT Growth relative to the local 
%                         classifications used cut-offs as described in supp 
%                         see Figure S7A of Kwon et al. Dev Cell 2015
%                         
%
% dispInRegion: (PARAM) : numeric 
%                         Default : 0.7 ( in um) 
%                         MT Cort Dist < which track will be automatically be classified
%                         as end-on. 
%                         see Figure S7A of Kwon et al. Dev Cell 2015
%
% dispDiscard : (PARAM) : numeric 
%                         Default : 0.3 (im um)
%                         MT Cort Dist < which the track will be discarded.                         
%                         see Figure S7A of Kwon et al. Dev Cell 2015
% 
% OUTPUT NAMING : 
% outputDirectory : (PARAM) : character  
%                          Default : Current Directory (pwd)
%                          path to output file minus the filename
% 
%
% 
% 
% 
% OUTPUT: 
%       toPlot  : struct with fields
%                 .info.names:    nGroupx1 cell 
%                                such that toPlot.names{iGroup} = 
%                                char array of that groups name  
%
%                 .info.projList : nGroupx1 cell
%                                 such that toPlot.projList{iGroup}{iProj} 
%                                 provides a char array of the project path
%                                 (subRegional in the mitotic case)
%                 .info.colors : nGroupx1 cell
%
%                 .info.grouping : a nProject x 1 double 
%                                 serving as a numeric grouping variable
%                                 for the projects. 
% 
%                 .(ip.Results.dataCollect) : an nGroupx1 cell
%                                  such that
%                                  toPlot.(ip.Results.dataCollect){iGroup}
%                                  provides a rxc double array of
%                                  measurements corresponding to the input
%                                  of dataCollect (either dwellTime or
%                                  percentLateral) where r is equal to the
%                                  maximum number of MT measurements per project  
%                                  per group and c is the number of
%                                  projects in the group. NaN is used as filler to 
%                                  per column if the number of MT per a given project 
%                                  is less than r. Designed for matlab
%                                  boxplot input. 
%                
%                                             
% 
%% Check input
ip = inputParser;

ip.CaseSensitive = false;

% type of meas
ip.addParameter('measToCollect','percentLateral'); % 'dwellTime' other option

% Parameters for the classification cut-offs
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x));
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x));
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x));
ip.addParameter('orientDiscard',150,@(x) isnumeric(x)); 

ip.addParameter('groupColors',[]);

ip.addParameter('outputDirectory',pwd);
ip.addParameter('filename',[]); 

ip.parse(varargin{:});
%% Set up
strMinDisp = num2str(ip.Results.minDispVect);
strOrient = num2str(ip.Results.orientCutOff);
strOrient = strrep(strOrient,'.','pt'); 
strDisp = num2str(ip.Results.dispInRegion);
strDisp = strrep(strDisp,'.','pt');
strDiscard = num2str(ip.Results.dispDiscard);
strDiscard = strrep(strDiscard,'.','pt');
strOrientDiscard = num2str(ip.Results.orientDiscard); 
fieldname = ['params' strMinDisp '_' strOrient '_' strDisp '_' strDiscard '_' strOrientDiscard];

if isempty(ip.Results.filename)
    filename = ['groupData_' fieldname '_' ip.Results.measToCollect];
else
    filename = ip.Results.filename;
end
names = unique(groupList(:,1),'stable');

if isempty(ip.Results.groupColors)
    nGroups = numel(names);
    colors = lines(nGroups);
    colors = arrayfun(@(x) colors(x,:) ,1:nGroups,'uniformoutput',0);
else
    colors = ip.Results.groupColors;
end

grpVar = cell(numel(names),1);
projList = cell(numel(names),1);
removeNoCortical = cell(numel(names),1); 
noViableData = cell(numel(names),1); 

%%
for iGroup = 1:numel(names)
    %Indices for current group
    idxC = cellfun(@(x) strcmpi(names{iGroup},x),groupList(:,1));
    projListC = groupList(idxC,2);
    
    idx =  arrayfun(@(x) exist([projListC{x} filesep 'meta' filesep 'CorticalInfo' filesep  'corticalData.mat'])==0,1:size(projListC,1));
    
    if sum(idx)>0
        idxRemove=  find(idx);
        arrayfun(@(x)  display(['No Cortical Info Exists for ' projListC{idxRemove(x)} ]),1:length(idxRemove));
        removeNoCortical{iGroup,1} = arrayfun(@(x) projListC{idxRemove(x)},1:length(idxRemove),'uniformoutput',0);
    else
        removeNoCortical{iGroup,1} = [];
    end
    
    % remove projects without cortical calculations. 
    projListC = projListC(~idx,:);
    
    % initiate dataMatC
    dataMatC = cell(length(projListC),1);
    
    
    for iProj = 1:numel(projListC)
        
            load([projListC{iProj} filesep 'meta' filesep 'projData.mat']);
            load([projListC{iProj} filesep 'meta' filesep 'CorticalInfo' filesep 'corticalData.mat']);
            
            switch ip.Results.measToCollect
                
                case 'dwellTime'
                    dataMatC{iProj,1} = corticalData.(fieldname).trackInfoEndOnST.dwells;
                    
                case 'percentLateral'
                    dataMatC{iProj,1} = corticalData.(fieldname).stats.percentLateral;
                    
                otherwise
                    error('Collection of that calculation is currently unsupported : check input or modify code');
            end
              
    end % for iProj
    
    
    % sometime have NO tracks that fall within the dwell category specified- can
    % either filter before or can catch here. - currently these will be
    % empty fields
    emptyTest  = cellfun(@(x) isempty(x),dataMatC);
    
    if sum(emptyTest) ~= 0
        idx = find(emptyTest);
        arrayfun(@(x) display(['No viable data from ' projListC{idx(x)} ': Removing']), 1:length(idx));
        noViableData{iGroup} = arrayfun(@(x) projListC{idx(x)},1:length(idx),'uniformoutput',0);
        
        % remove empty
        dataMatC = dataMatC(~emptyTest);
        projListC = projListC(~emptyTest);
    end
    
    projList{iGroup} = projListC;
    dataMatFinal = reformatDataCell(dataMatC);
    
    toPlot.(ip.Results.measToCollect).dataMat{iGroup} = dataMatFinal;
    grpVar{iGroup} = repmat(iGroup,size(projListC,1),1);  % 
end % iGroup

toPlot.info.grouping = vertcat(grpVar{:});
toPlot.info.colors = colors;
toPlot.info.projList = projList;
toPlot.info.names = names;

save([ip.Results.outputDirectory filesep filename ],'toPlot','removeNoCortical','noViableData'); 

end

