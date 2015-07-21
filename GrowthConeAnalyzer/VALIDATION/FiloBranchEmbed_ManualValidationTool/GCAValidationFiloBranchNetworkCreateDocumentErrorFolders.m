function [ output_args ] = GCAValidationReconstructMakeOverlays(selectedProjects,outDir)
%GCAValidationReconstructMakeOverlays

for iProj = 1:numel(selectedProjects)
    
    currentProj = selectedProjects{iProj,1};
    
    %get the ID  (find the function for this)
    [group,numID,date] = helperGetIDInfo(currentProj);
    outDirC = [outDir filesep num2str(iProj,'%03d')... 
        '_' date '_' group '_' numID '_Frame_' num2str(selectedProjects{iProj,2}) ];
    if ~isdir(outDirC)
        mkdir(outDirC);
    end 
    % load the MD file
    load([selectedProjects{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
    % make the troubleshoot reconstruction directory 
    GCATroubleshootMakeMovieOfReconstructMovie(MD,'frames',selectedProjects{iProj,2},...
        'outputDirectory',outDirC);
    
end 

end

