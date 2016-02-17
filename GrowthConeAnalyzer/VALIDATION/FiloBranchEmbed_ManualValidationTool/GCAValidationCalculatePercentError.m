function [ output_args ] = GCAValidationCalculatePercentError(reconstructDir,varargin)
%GCAValidationCalculatePercentError 

%% 
ip = inputParser;

ip.CaseSensitive = false;
ip.addRequired('reconstructDir');
ip.addParameter('type','filoBranch');

ip.parse(reconstructDir,varargin{:});

type = ip.Results.type;

x = dir(reconstructDir);
x = x(3:end);
list = arrayfun(@(i) x(i).name,1:length(x),'uniformoutput',0);

% will eventually have a filter here that will check if done and only load
% those not finished if user desires.


idxInclude  = listSelectGUI(list,[],'move');


% listSelectGUI to select a file to work on
toTest = list(idxInclude);

for iProj = 1: numel(toTest)
   
    switch type
        case 'filoBranch'
            load([reconstructDir filesep toTest{iProj} filesep ...
                'OverlayFiloBranch' filesep 'Reconstruct_Movie' filesep ...
                'filoInfo.mat']);
            FN = load([reconstructDir filesep toTest{iProj} filesep ...
                'OverlayFiloBranch' filesep 'ValidationOverlays' filesep ...
                'False_Neg_Coords.mat']);
            FP = load([reconstructDir filesep toTest{iProj} filesep ...
                'OverlayFiloBranch' filesep 'ValidationOverlays' filesep ...
                'False_Pos_Coords.mat']);
            add = []; 
            % currently for filoBranch I am saving a lot of information the
            %% if collect False Negatives 
            if ~isnan(FN.coords(1));
                
                for i  = 1:size(FN.coords,1)
                    x = arrayfun(@(x) ~isempty(intersect(sub2ind(imgSize,FN.coords(i,2),FN.coords(i,1)),filoInfo(x).Ext_maskIndices(:))),1:length(filoInfo));
                    
                end
            end
            
        case 'embed'
            load([reconstructDir filesep toTest{iProj} filesep ...
                'OverlayEmbed'  filesep 'filoInfo.mat']);
            FN = load([reconstructDir filesep toTest{iProj} filesep ...
                'OverlayEmbed' filesep 'ValidationOverlays' filesep ...
                'False_Neg_Coords.mat']);
            FP = load([reconstructDir filesep toTest{iProj} filesep ...
                'OverlayEmbed' filesep 'ValidationOverlays' filesep ...
                'False_Pos_Coords.mat']);
            add = 'Embedded'; 
            
            %MARIA CHECK THIS 
%              toRemove1 = arrayfun(@(x) numel(x.('Int_coordsXY')),filoInfoG);
%     toRemove2= arrayfun(@(x) isempty(x.('Int_coordsXY')),filoInfoG);
%     toKeep = ~(toRemove1==1 | toRemove2 ==1 );
%     NTot = sum(toKeep); 
    end
    
    % calculate value 
    NTot = length(filoInfoG); % filoInfoG should be the filtered value 
   
    %NTot = sum(toKeep); 
    
    nTrue(iProj,1) = NTot -size(FP.coords,1) + size(FN.coords,1);
    nPos(iProj,1) = size(FP.coords,1); 
    nNeg(iProj,1) = size(FN.coords,1); 
    percentPositive(iProj,1) = (size(FP.coords,1)/nTrue(iProj,1))*100; 
    percentNegative(iProj,1) = (size(FN.coords,1)/nTrue(iProj,1))*100; 
    
    
end
% rowNames = toTest; 
test = [percentPositive percentNegative nTrue nPos nNeg ]; 


varNames{1,1} = ['Percent_False_Positive_' add]; 
varNames{1,2} = ['Percent_False_Negative_' add]; 
varNames{1,3} = ['Ground_Truth_' add 'Filopodia'];
varNames{1,4} = ['Number_False_Positive_Filopodia' add]; 
varNames{1,5} = ['Number_False_Negative_Filopodia' add]; 

x = dataset({test,varNames{:}},'ObsNames',toTest); 
export(x,'File',['PercentError' add '.csv']); 
 
% False_PositiveEmbedded_Percent = percentPositive.*100; 
% False_NegativeEmbedded_Percent = percentNegative.*100;  

% t  = table([False_PositiveEmbedded_Percent,False_NegativeEmbedded_Percent],'RowNames',rowNames); 
% writetable(t,'PercentErrorEmbedded.csv'); 

