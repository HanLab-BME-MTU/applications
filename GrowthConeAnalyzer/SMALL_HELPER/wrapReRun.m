function [ output_args ] = wrapReRun(projList)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for iProj =1:numel(projList)
    %MD = load([projList{iProj} filesep 'movieData.mat']);
    
    to = [projList{iProj} filesep 'ANALYSIS_BEFORE_20140923'];  
    from = [projList{iProj} filesep 'ANALYSIS']; 
    
    movefile(from,to)
    
    % reRun Reconstruct
    cd([projList{iProj}]); 
    
 % performNeuriteOrientation   
    makeMovieDataMultChannels_NoBio(pwd,1,'neuriteOrient') 
    
    
    
  
    
end 
   

end

