function [ wekaDataset ] = getWekaDatasetFromCSVFile( csvFile )

    import weka.core.converters.*;
    import weka.core.*;
    import java.io.*;
    
    csvLoader = CSVLoader();
    csvLoader.setSource( File(java.lang.String(csvFile)) );
    wekaDataset = csvLoader.getDataSet();

end