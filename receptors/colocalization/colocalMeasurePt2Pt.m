function [cT, cNull, cCritical,estimatorM, estimatorC] = colocalMeasurePt2Pt(detectionRef,detectionObs, threshold, maskingFile)
%COLOCALMEASUREPT2PT measures the colocalization between two punctate channels 
%
%   Synopsis: [cT, cNull, cCritical,estimatorM, estimatorC] = colocalMeasurePt2Pt(detectionRef,detectionObs, threshold, maskingFile)
%
%   Function analyzes the distances of nearest neighbors between two
%   punctate channels and uses a distance threshold to ouput the probability of
%   finding points in the observed channel that colocalize with points in the
%   reference channel, the associated interaction potential (see source below)
%   and the null hypothesis values for both of these.
%Input:
% detectionRef: Detected positions of molecules in the reference channel
% 
%
% detectionObs: Detected positions of molecules in the observed channel
% 
%
% Threshold: distance threshold for possible interacting objects
%
% Masking File: This should either be a binary mask, or if no mask is given then the size of the image should be input
%               to create a mask of all ones
%
%
% Output:
% cNull (Null Hypothesis): probabilty of finding a distance <= threshold to
%   a detection in detectionRef under the null hypothesis (No interaction)
% %
% cT (Colocalization Measure): probabilty of finding a detection in detectionObs
%   a distance <= threshold to a detection in detectionRef
%
% estimatorC (Critical Potential): Minimum interaction potential require to reject
%   null hypothesis
%
% estimatorM (Calculated Potential): Calculated interaction potential from data 
%
% Estimator methodology adapted from "Beyond co-localization: inferring
% spatial interactions between sub-cellular structures from microscopy
% images" - Helmuth et al. BMC Bioinformatics 2010 
% Anthony Vega 09/2014

% Analysis
      
    %% Index points from both images
    % Populate query points and DT points with positions taken from images
    % Ref Channel
    xIndex = detectionRef.xCoord(:,1);
    yIndex = detectionRef.yCoord(:,1);
    QP = [yIndex xIndex];
    
    %Find detections in QP that lie inside maskList
    if ~ismatrix(maskingFile)
        maskingFile= ones(maskingFile);
    end
    [row, col] = find(maskingFile); 
    maskList = [row, col];
    lia1 = ismember(round(QP),maskList,'rows');
    
    %Multipling lia (binary vector) by QP will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    QP(:,1) = lia1.*QP(:,1);
    QP(:,2) = lia1.*QP(:,2);
    QP( ~any(QP,2), : ) = [];
    
% Obs Channel
    xIndex = detectionObs.xCoord(:,1);
    yIndex = detectionObs.yCoord(:,1);
    DT = [yIndex xIndex];
    
    %Find detections in QP that lie inside maskList
    [row, col] = find(maskingFile); 
    maskList = [row, col];
    lia1 = ismember(round(DT),maskList,'rows');
    
    %Multipling lia (binary vector) by QP will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    DT(:,1) = lia1.*DT(:,1);
    DT(:,2) = lia1.*DT(:,2);
    DT( ~any(DT,2), : ) = [];
    
    %Get nearest neighbor for every point in observed channel DT
    [~, d]= knnsearch(QP,DT,'K',1);
    
    % Test NN under threshold
%     xi = min(d):0.01:max(d); %possibly remove for consistency, maybe linspace
    [f,xi] = ksdensity(d);
    %Determine density of distances below threshold to find cNull
    test = find(xi<=threshold);
    cT = trapz(xi(1:max(test)),f(1:max(test))); 
    

    %% Test against null hypothesis
    [~, d1] = knnsearch(QP,maskList,'K',1);
 
    %Produce probablilty density of values
%     xi = min(d1):0.01:max(d1);
    
    [f,xi] = ksdensity(d1);
    %Determine density of distances below threshold to find cNull
    test = find(xi<=threshold);
    cNull(1,1) = trapz(xi(1:max(test)),f(1:max(test)));


    %Calculate estimator
    estimatorM = log(cT(1,1)/(1-cT(1,1)))-log(cNull(1,1)/(1-cNull(1,1)));


    % Estimate critical parameters based on cNull and interaction size
    cCritical = (binoinv(0.95,length(d),cNull(1,1)))/length(d);
    estimatorC = log(cCritical(1,1)/(1-cCritical(1,1)))-log(cNull(1,1)/(1-cNull(1,1)));
     

 end