%% Input: file name of the format: date (yymmdd), tumor # (m###m###), 'ratcolxy', exp # (nn)
% Input: date (6), tumor1-10(4), tumor11-20(4), **, exp # (2)
% Output: date, tumor (tumor string), tumorInd, isHighMetPot, day
function [expDetails] = parseExpName(expName,metaData)
expDetails.date = expName(1:6);

tumor1 = lower(expName(7:10));
tumor2 = lower(expName(11:14));

expDetails.expNum = str2num(expName(end-1:end));

if isempty(expDetails.expNum)
    return;
end


% Not in metaData or in corresponding exclude list -- exit
[tmp,indInMetaData] = ismember(expName(1:end-4),metaData.experiments.fnames);
if indInMetaData == 0 || ismember(expDetails.expNum,metaData.experiments.exclude{indInMetaData})
    return; 
end

% experiments 1-10 - tumor #1; experiments 11-20 - tumor #2;
if expDetails.expNum > 0 && expDetails.expNum <= 10
    expDetails.tumor = tumor1;
    %     expDetails.jobTumorPair = [num2str(indInMetaData) '_' tumor1];
    expDetails.jobTumorPair = [expDetails.date ' ' tumor1];    
else if expDetails.expNum > 10 && expDetails.expNum <= 20
        expDetails.tumor = tumor2;
%         expDetails.jobTumorPair = [num2str(indInMetaData) '_' tumor2];        
        expDetails.jobTumorPair = [expDetails.date ' ' tumor2];        
    else error('ilegel exp number\n');
    end
end

expDetails.jobTumorPairInd = find(ismember(metaData.experiments.jobTumorPair,expDetails.jobTumorPair));

expDetails.tumorInd = find(ismember(metaData.tumors.ids,expDetails.tumor));
expDetails.isHighMetPot = metaData.tumors.metastaticEfficiency(expDetails.tumorInd);

% TODO: days not implemented



%% OLD META DATA
% expDetails.date = expName(1:6);
% 
% tumor1 = expName(7:10);
% tumor2 = expName(11:14);
% 
% expDetails.expNum = str2num(expName(end-1:end));
% 
% if expDetails.expNum > 0 && expDetails.expNum <= 10
%     expDetails.tumor = tumor1;
% else if expDetails.expNum > 10 && expDetails.expNum <= 20
%         expDetails.tumor = tumor2;
%     else error('ilegel exp number\n');
%     end
% end
% 
% expDetails.tumorInd = find(ismember(metaData.tumors,expDetails.tumor));
% expDetails.isHighMetPot = find(ismember(metaData.highMetPot,expDetails.tumor));
% 
% if find(ismember(metaData.day1,expDetails.date))
%     expDetails.day = 1;
% else if find(ismember(metaData.day2,expDetails.date))
%         expDetails.day = 2;
%     else if find(ismember(metaData.day3,expDetails.date))
%             expDetails.day = 3;
%         else error('ilegal day number\n');
%         end
%     end
% end

end