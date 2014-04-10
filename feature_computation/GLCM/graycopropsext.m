function [ glcm_features ] = graycopropsext(glcm, varargin)

    % List of all GLCM properties this function can compute
    % Note for developers - if you want to implement a new 
    % property add a function to this file named ComputeGLCM_<property-name>
    % which computes the property and add <property-name> to this list
    glcmPropertyListSupported = { 'Contrast', ...
                                  'Dissimilarity', ... 
                                  'Homogeneity' ...
                                  'Energy', ...
                                  'Entropy', ...
                                  'Correlation', ...
                                  'AutoCorrelation', ...
                                  'ClusterShade', ...
                                  'ClusterProminence', ...                                  
                                  'MaximumProbability', ... 
                                  'SumOfSquaresVariance', ...
                                  'SumAverage', ...
                                  'SumVariance', ...
                                  'SumEntropy', ...
                                  'DifferenceEntropy', ...
                                  'DifferenceVariance' };
                        
    if nargin == 0
        glcmPropertyListSupported
        return;
    end

    % validate input arguments
    p = inputParser;
    p.addRequired('glcm', @(x) (isnumeric(x) && ndims(x) >= 2 && all(x(:)>=0)));
    p.addOptional('properties', 'basic', ...
                  @(x) ((ischar(x) && ismember(x, {'basic', 'all'})) || ...
                        (iscell(x) && all(ismember(x,glcmPropertyListSupported))) ));
    p.parse(glcm, varargin{:});
    
    % prepare the list of glcm properties that need to be computed
    glcmProperties = p.Results.properties;    
    if ischar(glcmProperties)
    
        switch glcmProperties
            
            case 'all'
                
                glcmProperties = glcmPropertyListSupported;
                
            case 'basic'
                
                glcmProperties = {'Contrast', ...
                                  'Homogeneity' ...
                                  'Energy', ...
                                  'Entropy', ...
                                  'Correlation', ...
                                  'AutoCorrelation', ...
                                  'ClusterShade', ...
                                  'ClusterProminence', ...
                                  'MaximumProbability' }; 
        end
        
    end
    
    numGLCM = size(glcm, 3);
    
    for g = 1:numGLCM 
        
        % get current glcm
        curGLCM = glcm(:, :, g);
        
        % normalize glcm
        normalized_glcm = curGLCM;

        if any(glcm(:))

            normalized_glcm = curGLCM ./ sum(curGLCM(:));

        end

        % Get row and column subscripts of GLCM.  These subscripts correspond 
        % to the pixel values in the GLCM.
        glcm_size = size(normalized_glcm);
        [c,r] = meshgrid(1:glcm_size(1), 1:glcm_size(2));
        r = r(:);
        c = c(:);

        % compute features
        for pid = 1:numel(glcmProperties)           
            propfunc = str2func( ['ComputeGLCM_' glcmProperties{pid}] );
            cur_glcm_features.(glcmProperties{pid}) = propfunc(normalized_glcm, r, c);            
        end

        % add current glcm to list
        glcm_features(g) = cur_glcm_features;
        
    end
            
end

%% GLCM Contrast Feature Group

% GLCM Contrast ---- sum(Cij x (i - j)^2)
function [ Contrast ] = ComputeGLCM_Contrast(normalized_glcm, r, c)

    term1 = (r - c).^2;
    term2 = normalized_glcm(:);

    Contrast = sum(term1 .* term2);

end

% GLCM Dissimilarity
function [ Dissimilarity ] = ComputeGLCM_Dissimilarity(normalized_glcm, r, c)

    term1 = abs(r - c);
    term2 = normalized_glcm(:);

    Dissimilarity = sum(term1 .* term2);

end

% GLCM Homogeneity
function [ Homogeneity ] = ComputeGLCM_Homogeneity(normalized_glcm, r, c)

    term1 = 1 + (r - c).^2;
    term2 = normalized_glcm(:);

    Homogeneity = sum(term2 ./ term1);
    
end

%% GLCM Orderliness Feature Group

% GLCM Angular Second Moment (ASM) : sum(Cij.^2)
function [ ASM ] = ComputeGLCM_AngularSecondMoment(normalized_glcm, r, c)

    term1 = normalized_glcm(:);

    ASM = sum(term1.^2);
    
end

% GLCM Energy : sqrt(ASM)
function [ Energy ] = ComputeGLCM_Energy(normalized_glcm, r, c)

    ASM = ComputeGLCM_AngularSecondMoment(normalized_glcm, r, c);
    
    Energy = sqrt(ASM);
    
end

% GLCM Entropy : -1 x sum(Cij x log(Cij))
function [ Entropy ] = ComputeGLCM_Entropy(normalized_glcm, r, c)

    term1 = normalized_glcm(:);    
    
    term2 = log(term1 + eps);
    
    Entropy = -1.0 * sum(term1 .* term2);
    
end

%% GLCM Stats Group

% GLCM Marginal Mean
function [ M ] = ComputeGLCM_Mean_Index(normalized_glcm, index)

    M = sum(index .* normalized_glcm(:));

end

% GLCM Marginal Variance
function [ V ] = ComputeGLCM_Variance_Index(normalized_glcm, index)

    M = sum(index .* normalized_glcm(:));
    
    V = sum(((index - M).^2) .* normalized_glcm(:));
    
end

% GLCM Auto Correlation
function [ corr ] = ComputeGLCM_AutoCorrelation(normalized_glcm, r, c)

    corr = sum(normalized_glcm(:) .* r .* c);

end

% GLCM Correlation
function [ C ] = ComputeGLCM_Correlation(normalized_glcm, r, c)

    mr = ComputeGLCM_Mean_Index(normalized_glcm, r);
    mc = ComputeGLCM_Mean_Index(normalized_glcm, c);
    
    sr = sqrt(ComputeGLCM_Variance_Index(normalized_glcm, r));
    sc = sqrt(ComputeGLCM_Variance_Index(normalized_glcm, c));
    
    term1 = ((r - mr) .* (c - mc));
    term2 = sr * sc;
    
    if term2 == 0        
        
        C = 1;    % variance is zero        
        
    else
        
        C = sum(normalized_glcm(:) .* (term1 ./ term2));        
        
    end    
    
end

% GLCM Cluster Shade : sum(Cij * ((i - mi) + (j - mj)).^3)
function [ Shade ] = ComputeGLCM_ClusterShade(normalized_glcm, r, c)

    mr = ComputeGLCM_Mean_Index(normalized_glcm, r);
    mc = ComputeGLCM_Mean_Index(normalized_glcm, c);    
    
    term1 = ((r - mr) + (c - mc)).^3;
    term2 = normalized_glcm(:);
    
    Shade = sum(term1 .* term2);
    
end

% GLCM Cluster Prominence : sum(Cij * ((i - mi) + (j - mj)).^4)
function [ Prominence ] = ComputeGLCM_ClusterProminence(normalized_glcm, r, c)

    mr = ComputeGLCM_Mean_Index(normalized_glcm, r);
    mc = ComputeGLCM_Mean_Index(normalized_glcm, c);    
    
    term1 = ((r - mr) + (c - mc)).^4;
    term2 = normalized_glcm(:);
    
    Prominence = sum(term1 .* term2);
    
end

%% Miscelleneous (all other types of) features group
function [ maxProb ] = ComputeGLCM_MaximumProbability(normalized_glcm, r, c)

    maxProb = max(normalized_glcm(:));

end

% GLCM Sum of Squares Variance (to me this seems like an absurd feature
% r is gray-level and mu is probability, why wud any sane soul take a 
% design a feature as function of the difference between them)
function [ svar ] = ComputeGLCM_SumOfSquaresVariance(normalized_glcm, r, c)

    mu = mean2(normalized_glcm);    
    svar = sum(normalized_glcm(:) .* (r - mu).^2);
    
end

% GLCM Sum Average
function [ savg ] = ComputeGLCM_SumAverage(normalized_glcm, r, c)

    k = r + c - 2;
    savg = sum(k .* normalized_glcm(:));
    
end

% GLCM Sum Entropy
function [ sument ] = ComputeGLCM_SumEntropy(normalized_glcm, r, c)

    numGrayLevels = size(normalized_glcm,1);
    
    indsum = r + c - 2;   
    
    psum = zeros(2*numGrayLevels-1,1);
    for k = 0:(2*(numGrayLevels-1))
        psum(k+1) = sum(normalized_glcm(indsum == k));
    end
    
    sument = -1.0 * sum( psum .* log(psum + eps) );
    
end

% GLCM Sum Variance
function [ sumvar ] = ComputeGLCM_SumVariance(normalized_glcm, r, c)

    numGrayLevels = size(normalized_glcm,1);
    
    indsum = r + c - 2;   
    
    psum = zeros(2*numGrayLevels-1,1);
    kvals = 0:(2*(numGrayLevels-1));
    for k = kvals
        psum(k+1) = sum(normalized_glcm(indsum == k));
    end
    
    sument = -1.0 * sum( psum .* log(psum + eps) );    
    sumvar = sum((kvals' - sument).^2 .* psum);
    
end

% GLCM Difference Entropy
function [ diffent ] = ComputeGLCM_DifferenceEntropy(normalized_glcm, r, c)

    numGrayLevels = size(normalized_glcm,1);
    
    inddiff = abs(r - c);   
    
    pdiff = zeros(numGrayLevels,1);
    for k = 0:(numGrayLevels-1)
        pdiff(k+1) = sum(normalized_glcm(inddiff == k));
    end
    
    diffent = -1.0 * sum( pdiff .* log(pdiff + eps) );
    
end

% GLCM Difference Variance
function [ diffvar ] = ComputeGLCM_DifferenceVariance(normalized_glcm, r, c)

    numGrayLevels = size(normalized_glcm,1);
    
    inddiff = abs(r - c);   
    
    pdiff = zeros(numGrayLevels,1);
    kvals = 0:(numGrayLevels-1);
    for k = kvals
        pdiff(k+1) = sum(normalized_glcm(inddiff == k));
    end
    
    diffvar = sum( (kvals').^2 .* pdiff );
end

