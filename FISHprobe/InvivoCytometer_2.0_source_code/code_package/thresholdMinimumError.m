function [thresholdValue, varargout] = thresholdMinimumError(imInput, varargin)
% % Computes a threshold using the minimum error thresholding method (See Ref [1]).
% 
%     [ thresholdValue ] = thresholdMinimumError(imInput, varargin );
%     [ thresholdValue, imMask] = thresholdMinimumError(imInput,varargin );
%     [ thresholdValue, imMask, minError ] = thresholdMinimumError(imInput,varargin );
%     
% The image histogram is modeled as a mixture of two gaussian/poisson 
% distributions - one for the background and one for the foreground. 
% The threshold is then computed by minimizing the relative entropy
% between histogram and the mixture model. 
% 
% It is common/popular to use a gaussian model for the mixture components, 
% but there is a line of thought that advocates the use of a poisson model.
% 
% This function assumes a bright foreground against a relatively darker 
% background.
% 
%     Required Input Arguments:     
%         
%                          imInput: Input ND Image
%                          
%     Optional Input Arguments:    
% 
%                 numHistogramBins: number of bins to use for the histogram                
%                                   Default: 256
%                                   
%     Optional Param/Value Arguments:    
%                                
%                            model: specifies the distribution should be used 
%                                   to model the foreground and background 
%                                   portion of the histogram.
%                                   Options: 'gaussian' (popular choice)
%                                            'poisson'  
%                                   Default: 'poisson'
%                                   
%                        debugMode: true/false
%                                   specifies whether or not to run in 
%                                   debug mode. A bunch of stuff is 
%                                   printed/plotted. This mode is only 
%                                   intended for the developers of this 
%                                   function. If you are not one, dont
%                                   bother.                                  
%                                   Default: false
%                                   
% Examples:
% 
%     I = imread('coins.png');
%     
%     figure;
%     
%     subplot(2,2,1)
%     
%         imshow( I );
%         title( 'Input Image', 'FontWeight', 'bold' );
%         
%     subplot(2,2,2)
%     
%         levelOtsu = graythresh(I);
%         imshow( im2bw(I, levelOtsu) );
%         title( sprintf('Otsu Threshold - %d', round(levelOtsu*255) ), ...
%                'FontWeight', 'bold' );
% 
%     subplot(2,2,3)
%     
%         levelMinErrPoisson = thresholdMinimumError( I, 'model', 'poisson' );
%         imshow( I >= levelMinErrPoisson );
%         title( sprintf('Poisson Min-error Threshold - %d', round(levelMinErrPoisson) ), ...
%                'FontWeight', 'bold' );
% 
%     subplot(2,2,4)
%     
%         levelMinErrGaussian = thresholdMinimumError( I, 'model', 'gaussian' );
%         imshow( I >= levelMinErrGaussian );
%         title( sprintf('Gaussian Min-error Threshold - %d', round(levelMinErrGaussian) ), ...
%                'FontWeight', 'bold' );
%         
% References:
% 
% [1] Fan, J. (1998). "Notes on Poisson distribution-based minimum error 
% thresholding." Pattern Recognition Letters 19(5–6): 425-431.
% 
% [2] PAL, N. R. and S. K. PAL (1991). "IMAGE MODEL, POISSON DISTRIBUTION 
% AND OBJECT EXTRACTION." International Journal of Pattern Recognition and
% Artificial Intelligence 05(03): 459-483.
% 
% Author: Deepak roy Chittajallu (Created Feb 28, 2013)
% 

    p = inputParser;
    p.addRequired('imInput', @isnumeric);
    p.addOptional('numHistogramBins', 256, @(x) (isscalar(x) && (x-floor(x)) == 0));
    p.addParamValue('model', 'gaussian', @(x) (ismember(x, {'gaussian', 'poisson'})));
    p.addParamValue('debugMode', false, @islogical);
    p.parse(imInput, varargin{:});
    
    numHistogramBins = p.Results.numHistogramBins;
    probModel = p.Results.model;
    flagDebugMode = p.Results.debugMode;

    imInputOrg = imInput;
    
    % Convert ND image to 1D (should makes things faster) and cast to double
    imInput = double(imInput(:));
    
    % compute image intensity range
    minIntensity = min(imInput);
    maxIntensity = max(imInput);
    
    % normalize image to [0,1]
    grayToHistMultiplier = numHistogramBins / (maxIntensity - minIntensity);
    imInput = (imInput - minIntensity) * grayToHistMultiplier;
    
    % compute normalized histogram
    [binFreq, grayVal] = hist(imInput, numHistogramBins);
    
    binFreq = binFreq / numel(imInput);
    binFreq = (binFreq + eps);
    binFreq = binFreq / sum(binFreq);
    
    binFreqCumSum = cumsum(binFreq);    
    trueGrayVal = minIntensity + grayVal / grayToHistMultiplier;
    
    % compute the threshold that minimizes the relative entropy between 
    % the model and the histogram
    errorFunction = zeros(size(binFreq));    
    gCumSum = cumsum(grayVal .* binFreq);
    
    switch probModel
        
        case 'poisson'            

            priorLeft = binFreqCumSum;
            priorRight = binFreqCumSum(end) - binFreqCumSum;
            
            meanWhole = gCumSum(end);
            meanLeft = gCumSum ./ priorLeft;
            meanRight = (gCumSum(end) - gCumSum) ./ priorRight;
            
            errorFunction = meanWhole ...
                            - priorLeft .*  (log(priorLeft) + (meanLeft .* log(meanLeft))) ...
                            - priorRight .* (log(priorRight) + (meanRight .* log(meanRight)));
                           
        case 'gaussian'            

            gSqCumSum = cumsum((grayVal.^2) .* binFreq);

            priorLeft = binFreqCumSum;
            priorRight = binFreqCumSum(end) - binFreqCumSum;

            meanWhole = gCumSum(end);
            meanLeft = gCumSum ./ priorLeft;
            meanRight = (meanWhole - gCumSum) ./ priorRight;

            % variance = E[X^2] - (E[X])^2
            varLeft = (gSqCumSum ./ priorLeft) - meanLeft.^2;
            varRight = ((gSqCumSum(end) - gSqCumSum) ./ priorRight) - meanRight.^2;

            stdLeft = sqrt(varLeft) + eps;
            stdRight = sqrt(varRight) + eps;

            errorFunction = ((priorLeft .* log(stdLeft)) + (priorRight .* log(stdRight))) ...
                            - ((priorLeft .* log(priorLeft)) + (priorRight .* log(priorRight)));
                        
            errorFunction([1,end]) = max(errorFunction(2:end-1)) + eps;
    end

    [minError, minInd] = min( errorFunction );
    thresholdValue = trueGrayVal(minInd);
    
    if flagDebugMode
        
        figure;
        
        subplot(2,1,1)
        
            plot(trueGrayVal, errorFunction, 'b-', 'LineWidth', 2.0);
            title( sprintf( 'Optimal Threshold = %f', thresholdValue ), ...
                   'FontWeight', 'bold' );
            xlabel( 'Intensity' );
            ylabel( 'Error' );
            
            hold on;
                plot( trueGrayVal(minInd), minError, 'mx', 'MarkerSize', 15.0 );
            hold off;
        
        subplot(2,1,2) 
        
            bar(trueGrayVal, binFreq);
            ybounds = ylim;
            hold on;
            
                % show optimal threshold
                plot( trueGrayVal(minInd)*ones(1,2), [0,ybounds(2)], 'm-', 'LineWidth', 2.0 );
                
                % overlay the foreground background pdfs
                switch probModel
                
                    case 'poisson'
                        
                        % estimate parameters
                        priorLeft = binFreqCumSum(minInd);
                        priorRight = binFreqCumSum(end) - binFreqCumSum(minInd);

                        meanWhole = gCumSum(end);
                        meanLeft = (gCumSum(minInd) / priorLeft);
                        meanRight = ((meanWhole - gCumSum(minInd)) / priorRight);
                        
                        % compute probabilities for whole intensity range
                        bgProb = poisspdf(1:numHistogramBins, meanLeft);
                        bgProb = bgProb / sum(bgProb);

                        fgProb = poisspdf(1:numHistogramBins, meanRight);
                        fgProb = fgProb / sum(fgProb);
                        
                        % plot pdfs
                        plot( trueGrayVal, bgProb, 'r-', 'LineWidth', 2.0 );
                        plot( trueGrayVal, fgProb, 'g-', 'LineWidth', 2.0 );
                        
                    case 'gaussian'    
                   
                        % estimate parameters
                        priorLeft = binFreqCumSum(minInd);
                        priorRight = binFreqCumSum(end) - binFreqCumSum(minInd);

                        meanWhole = gCumSum(end);
                        meanLeft = (gCumSum(minInd) / priorLeft);
                        meanRight = ((meanWhole - gCumSum(minInd)) / priorRight);

                        varLeft = (gSqCumSum(minInd) / priorLeft) - meanLeft^2;
                        varRight = ((gSqCumSum(end) - gSqCumSum(minInd)) / priorRight) - meanRight^2;

                        stdLeft = sqrt(varLeft);
                        stdRight = sqrt(varRight);
                        
                        % compute probabilities for whole intensity range
                        bgProb = normpdf(1:numHistogramBins, meanLeft, stdLeft);
                        bgProb = bgProb / sum(bgProb);

                        fgProb = normpdf(1:numHistogramBins, meanRight, stdRight);
                        fgProb = fgProb / sum(fgProb);
                        
                        % plot pdfs
                        plot( trueGrayVal, bgProb, 'r-', 'LineWidth', 2.0 );
                        plot( trueGrayVal, fgProb, 'g-', 'LineWidth', 2.0 );
                        
                end
                
            hold off;
        
    end
    
    if nargout > 1
        varargout{1} = imInputOrg > thresholdValue;        
    end    
    
    if nargout > 2
       varargout{2} = minError;
    end
    
end