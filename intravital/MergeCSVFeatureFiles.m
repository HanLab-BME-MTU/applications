

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    classificationTask = { 'G1_S_G2M', 'G2_M' };

    defaultOutputDir = 'Z:\intravital\data\Stefan_September2012 (Mouse 4)\imageanalysis';    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ask the user to provide the list of directories to be processed
rootDirList = uipickfiles('FilterSpec', defaultOutputDir);

% ask the user to select the output directory
outputRootDir = uigetdir( defaultOutputDir, 'Select Output Directory' ); 

% merge csv files pairwise
AddWekaClassesToPath();

import weka.core.*;

filecombs = nchoosek( rootDirList, 2 );

for i = 1:size(filecombs,1)
    
    fprintf('\n\nGenerating Merged Feature File - %d/%d\n\n', i, size(filecombs,1));
    
    [pathstr, dirName1, ext] = fileparts(filecombs{i,1});
    [pathstr, dirName2, ext] = fileparts(filecombs{i,2});
    
    for cid = 1:numel(classificationTask)

        % copy individual csv files to output directory
        featureFile1_src = fullfile( filecombs{i,1}, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
        featureFile1_dest = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s_%s.csv', classificationTask{cid}, dirName1));
        copyfile(featureFile1_src, featureFile1_dest, 'f');
        ConvertCSVFileToArffFile( featureFile1_dest );
        
        featureFile2_src = fullfile( filecombs{i,2}, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
        featureFile2_dest = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s_%s.csv', classificationTask{cid}, dirName2));
        copyfile(featureFile2_src, featureFile2_dest, 'f');        
        ConvertCSVFileToArffFile( featureFile2_dest );
        
        % append the two files and copy to output directory
        outFile = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s_%s_%s.csv', classificationTask{cid}, dirName1, dirName2));
        AppendTwoFeatureFiles( featureFile1_src, featureFile2_src, outFile );
        ConvertCSVFileToArffFile( outFile );
        
        
    end
    
end




