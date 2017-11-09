function [ success ] = GCAVisualsWrapper(projList,visFuncH,varargin)
%GCAVisualsWrapper: a wrapper function that takes in a list of projects and
% a visualization function handle and makes these graphics
%
% projList: cell array of projects
% visFuncH: (char) name of the visualization function to run

%% EXTRA check if run
% Already checked if run in previous step : add later if you want as a
% standalone part of the function
%idxToRun = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep  'masksOtsu' ])==0, projList);
%idxToRun{1} = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'Raw'])==0,projList);
%idxToRun{2} = cellfun(@(x) exist([x filesep 'ANALYSIS' filesep 'Visualization_Overlays' filesep 'Overlaysmono'
%projListRun = projList(~idxToRun{1});
%cellfun(@(x) display(['Project ' x ' Already Run: SKIPPING']),projListRun);
%projListToRun = projList(idxToRun);
%%
ip = inputParser; 
ip.CaseSensitive = false;
ip.addParameter('ScaleBar',true); 
ip.addParameter('Timer',true); 
ip.addParameter('TreatmentFrame',[]); 
ip.addParameter('TreatmentName','CK666'); 
ip.addParameter('Overwrite',false); 
ip.addParameter('NeuriteID',true); 


ip.parse(varargin{:});
p = ip.Results; 
%%

projListToRun = projList;

if ~isempty(projListToRun)
    for iProj = 1:size(projListToRun,1)
        display(['Making Visualization Overlays for ' projListToRun{iProj} 'using' visFuncH])  
        if ip.Results.NeuriteID 
           nameID =  strrep(projList{iProj,2},'_',' '); 
           nameID = strrep(nameID,'KDNo','Control'); 
           p.NeuriteIDText = nameID; 
        end 
        load([projListToRun{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
        % check if run
        handleC= str2func(visFuncH);
        
        success(iProj) = handleC(MD,p); % success will be 1 if the information is sufficient to run 
        % the visualization and 0 if not 
    end
else 
    success = NaN; % noProjects 
end
end

