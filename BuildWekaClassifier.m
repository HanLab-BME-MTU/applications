function [cls] = BuildWekaClassifier( wekaTrainingDataset, outFile, varargin )

    % parameters
    p = inputParser;
    
        % number of components in the ensemble
        p.addParamValue( 'numEnsembleComponents', 5, @(x) (isscalar(x) && isnumeric(x)) );

        % decision fusion rule for the ensemble of classifiers
        % must be a string with one of the following values:
        % AVG | PROD | MAJ | MIN | MAX | MED
        p.addParamValue( 'fusionRule', 'PROD', @(x) ( ischar(x) && ismember(x, {'AVG', 'PROD', 'MAJ', 'MIN', 'MAX', 'MED'}) ) );
        
        % number of trees used for random forest
        p.addParamValue( 'numTrees', 100, @(x) (isscalar(x) && isnumeric(x)) );

        % number of folds used for cross validation
        p.addParamValue( 'numCVFolds', 10, @(x) (isscalar(x) && isnumeric(x)) );
    
    p.parse( varargin{:} );
    PARAMETERS = p.Results;
    
    % weka imports
    import weka.core.*;
    import weka.core.converters.*;

    import weka.classifiers.meta.*;
    import weka.filters.supervised.instance.*;
    import weka.classifiers.functions.*;              
    import weka.classifiers.trees.*;       
    import weka.classifiers.*;
    import weka.classifiers.evaluation.output.prediction.*; 
    
    import java.util.*;
    import java.lang.*;

    % build the classifier     
    fprintf('\nGenerating the classifier ...\n' );
    
    cls = weka.classifiers.meta.Vote();

    strFusionRule = sprintf( '-R %s', PARAMETERS.fusionRule );
    cls.setOptions( java.lang.String(strFusionRule).split(' ') );

    ensembleComponents = [];
    for i = 1:PARAMETERS.numEnsembleComponents
        
        % create spread subsample filter
        classBalancingFilter = SpreadSubsample();                    
        classBalancingFilter.setDistributionSpread( 1.0 );
        classBalancingFilter.setRandomSeed( i );

        % create classifier
        clsComp = RandomForest();
        clsComp.setNumTrees( PARAMETERS.numTrees );

        % create filtered classfier
        balancedClassifier = FilteredClassifier();
        balancedClassifier.setFilter( classBalancingFilter );
        balancedClassifier.setClassifier( clsComp );

        % add classifier to ensemble
        ensembleComponents = [ensembleComponents; balancedClassifier];
        
    end

    cls.setClassifiers( ensembleComponents );
    cls.buildClassifier( wekaTrainingDataset );
    
    SerializationHelper.writeAll(outFile, cls);         

    % perform cross validation
    fprintf('\nPerforming cross validation ...\n' );
    
    clsOutput = PlainText();
    clsOutputBuf = java.lang.StringBuffer();
    clsOutput.setBuffer(clsOutputBuf);
    
    evaluator = weka.classifiers.Evaluation(wekaTrainingDataset);
    evaluator.crossValidateModel(cls, wekaTrainingDataset, PARAMETERS.numCVFolds, java.util.Random(1), clsOutput);
    
    [pathstr, name, ext] = fileparts( outFile );
    
    fidPerformance = fopen( fullfile(pathstr, 'model_performance.txt'), 'w' );    
    fidVec = [1, fidPerformance];
    
    fprintfc(fidVec, 'Classifier used:\n\n%s\n', char(cls.toString()) );
    fprintfc(fidVec, '\n%d-fold cross-validation is used to evaluate performance\n\n%s\n', PARAMETERS.numCVFolds);
    fprintfc(fidVec, '\nResults Summary:\n\n%s\n', char(evaluator.toSummaryString()) );
    fprintfc(fidVec, '\nDetailed Accuracy By Class:\n\n%s\n', char(evaluator.toClassDetailsString()) );
    fprintfc(fidVec, '\nConfusion Matrix:\n\n%s\n', char(evaluator.toMatrixString()) );
    
    GMOR = 1;
    avgFmeasure = 0;
    numClasses = wekaTrainingDataset.numClasses();
    
    for cid = 1:numClasses
        p = evaluator.precision(cid-1);
        r = evaluator.recall(cid-1);
        GMOR = GMOR * r;
        avgFmeasure = avgFmeasure + harmmean( [p, r] );
    end
    
    GMOR = GMOR^(1/numClasses);
    avgFmeasure = avgFmeasure / numClasses;
    
    fprintfc(fidVec, '\nGeometric Mean of Recalls - %f', GMOR );     
    fprintfc(fidVec, '\nAverage of Per-class F-Measure - %f\n', avgFmeasure );     
    
    fclose( fidPerformance );
    
end