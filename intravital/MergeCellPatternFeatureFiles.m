

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

for cid = 1:numel(classificationTask)

    fprintf('\n\nGenerating Merged Feature File for Classification Task - %s\n\n', classificationTask{i});
    
    featureFile_Out = fullfile(outputRootDir, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
    
    for i = 1:numel(rootDirList)

        % copy individual csv files to output directory
        if i == 1
            featureFile_src = fullfile( rootDirList{i}, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
            copyfile(featureFile_src, featureFile_Out, 'f');
        else
            featureFile_src = fullfile( rootDirList{i}, sprintf( 'cellCycleFeatures_%s.csv', classificationTask{cid}) );
            
            % append file
            fidOut = fopen(featureFile_Out, 'a');            
            fidAppend = fopen(featureFile_src, 'r'); 
            
            strLine = fgetl(fidAppend); % ignore first line with col headers

            strLine = fgetl(fidAppend);
            while ischar(strLine)
               fprintf(fidOut, '\n%s', strLine); 
               strLine = fgetl(fidAppend);
            end

            fclose(fidOut);            
            fclose( fidAppend );
            
        end
        
    end

    ConvertCSVFileToArffFile( featureFile_Out );
    
end




