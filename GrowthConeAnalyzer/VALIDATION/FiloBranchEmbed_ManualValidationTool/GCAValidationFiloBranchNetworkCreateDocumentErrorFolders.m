function [ output_args ] = GCAValidationReconstructMakeOverlays(selectedProjects,outDir)
%GCAValidationReconstructMakeOverlays

for iProj = 1:size(selectedProjects,1)
    
    currentProj = selectedProjects{iProj,1};
    
    %get the ID  (find the function for this)
    %[group,numID,date] = helperGetIDInfo(currentProj);
    
       
            [~,date] = upDirectory(currentProj,2,1);
            [~,numID] = upDirectory(currentProj,1,1);
            [~,group] = upDirectory(currentProj,3,1);
           
    
    outDirC = [outDir filesep num2str(iProj,'%03d')... 
        '_' date '_' group '_' numID '_Frame_' num2str(selectedProjects{iProj,2}) ];
    if ~isdir(outDirC)
        mkdir(outDirC);
    end 
    % load the MD file
    load([selectedProjects{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
    inDir = [MD.outputDirectory_ filesep '/SegmentationTesting/' ... 
        'filoBranchParamScan/scanParam_maxRadiusConnectFiloBranch/015/VII_filopodiaBranch_fits/Channel_1']; 
    
    
    
     
    % make the troubleshoot reconstruction directory 
    GCATroubleshootMakeMovieOfReconstructMovie(MD,'frames',selectedProjects{iProj,2},...
        'OutputDirectory',outDirC,'InputDirectory',inDir,'outDirType',[]);
    
end 

end

