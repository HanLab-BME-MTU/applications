function ConvertCSVFileToArffFile( csvFilePath, arffFilePath )

    if ~exist( 'arffFilePath', 'var' )
        [pathstr, name, ext] = fileparts( csvFilePath );
        arffFilePath = fullfile(pathstr, [name '.arff']);
    end
    
    wekaDataset = getWekaDatasetFromCSVFile( csvFilePath );
    writeWekaDatasetToArffFile( wekaDataset, arffFilePath );

end