function writeWekaDatasetToArffFile( wekaDataset, arffFilePath )

    import weka.core.converters.*;
    import weka.core.*;
    import java.io.*;
    
    arffFileSaver = ArffSaver();
    arffFileSaver.setFile( File(java.lang.String(arffFilePath)) );
    arffFileSaver.setInstances( wekaDataset );
    arffFileSaver.writeBatch();
    
end