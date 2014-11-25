function [ output ] = quantDwellStats( list  )

% make titles 
types{1} = 'Terminal';
types{2} = 'Pause';
types{3} = 'Shrinkage';
types{4} = 'Undefined';

paramNames{1} = 'Number Sampled';
paramNames{2} = 'Number With Co-Decay Criteria';
paramNames{3} =  'Percentage Co-Decay Criteria'; 
count =2 ;
for iType = 1:4
for iParam = 1:3 
output{count,1} = [paramNames{iParam} ' ' types{iType}]; 
count = 1+count; 
end 
end 
%% 
for iProj = 1:numel(list) 

    % load the dwellInt 
    load([list{iProj} filesep 'dwellInt.mat']); 
   [~,projName] = upDirectory(list{iProj},2,1);
   
   output{1,iProj+1} = projName; 

count = 2  ;
for iType = 1:4
    if ~isempty(dwellInt.consistencyTest{iType})
        
number = length(horzcat(dwellInt.consistencyTest{iType}{:})); 
 decay = sum(vertcat(dwellInt.consistencyTest{iType}{:})); 
percent = decay/number*100;
    else
        number = NaN;
        percent = NaN;
        decay= NaN; 
    end
    
    output(count:count+2,iProj+1) = [{number} ;  {decay} ; {percent}];
    count = count+3; 
end 
   
end 
    
   



end

 

