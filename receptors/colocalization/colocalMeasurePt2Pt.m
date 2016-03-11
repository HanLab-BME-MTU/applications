function [cT, cNull, cCritical,estimatorM, estimatorC] = colocalMeasurePt2Pt(detectionRef,detectionObs, threshold, maskingFile)
%COLOCALMEASUREPT2PT measures fraction of coloclization between two punctate channels 
% Measures the fraction of the detectionObs that colocalizes with
% detectionRef.
%Input:
% imageInfoA: Detected positions of molecules in channel which will be assigned a nearest neighbor
% 
%
% imageInfoB: Detected positions of molecules in second channel
% 
%
% Threshold: distance threshold for possible interacting objects
%
% Masking File: This should either be a binary mask, or is no mask is given then the size of the image should be input
%               to create a mask of all ones
%
%
% Output:
% cNull (Null Hypothesis): probabilty of finding a distance <= threshold among
% Nearest Neighbors (NN)
% %
% cT (Colocalization Measure): Fraction of NN objects which meet threshold
%
% estimatorC (Critical Potential): Minimum allowed potential given c0
%
% estimatorM (Calculated Potential): Calculated potential from c0 and c1


% Populate query points and DT points with positions taken from images
% xIndex: Information from detection program which reads x position and
% variance of objects in a single frame
%
% yIndex: Information from detection program which reads y position and
% variance of objects in a single frame
 a = 1;
      
    %Index points from both image
    xIndex = detectionRef.xCoord(:,1);
    yIndex = detectionRef.yCoord(:,1);
    QP = [yIndex xIndex];
    
    %Find detections in QP that lie inside maskList
    if ~ismatrix(maskingFile)
        maskingFile= ones(maskingFile);
    end
    [row, col] = find(maskingFile); %Verify
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
    [row, col] = find(maskingFile); %Verify
    maskList = [row, col];
    lia1 = ismember(round(DT),maskList,'rows');
    
    %Multipling lia (binary vector) by QP will replace coord outside
    %boundary with zero, last line removes all zeros from vector
    DT(:,1) = lia1.*DT(:,1);
    DT(:,2) = lia1.*DT(:,2);
    DT( ~any(DT,2), : ) = [];
    
    [~, d]= knnsearch(QP,DT,'K',1);
    
% Other stuff to be changed------------------------------------------------
    % Test NN under threshold
    xi = min(d):0.01:max(d);
    [f] = ksdensity(d,xi);
    %Determine density of distances below threshold to find cNull

    test = find(xi<=threshold);
    cT = trapz(xi(1:max(test)),f(1:max(test))); 
    
% %     figure; plot(xi,f);

    %Determine Ct score
    %cT(a,1) = length(C)/length(d);

    %Create another channel to measure q(d)
%     [~, d1] = nearestNeighbor(DT, sDensity);
%     A = [1:0.25:512];
%     Y = repmat(A,1,length(A));
%     X = repmat(A,length(A),1);
%     X = X(1:end);
%     X1(:,1) = X;
%     X1(:,2) = Y;
   [~, d1] = knnsearch(QP,maskList,'K',1);
 
    %Produce probablilty density of values
    xi = min(d1):0.01:max(d1);
    
    [f] = ksdensity(d1,xi);
    %Determine density of distances below threshold to find cNull
    test = find(xi<=threshold);
    cNull(a,1) = trapz(xi(1:max(test)),f(1:max(test)));

    % Plot pdf
    % figure()
    % plot(xi,f);
    % xlim([-5 50])
    % ylim([0 1])



    %Calculate estimator
    estimatorM = log(cT(a,1)/(1-cT(a,1)))-log(cNull(a,1)/(1-cNull(a,1)));


    % Estimate critical parameters based on cNull and interaction size
    cCritical = (binoinv(0.95,length(d),cNull(a,1)))/length(d);
    estimatorC = log(cCritical(a,1)/(1-cCritical(a,1)))-log(cNull(a,1)/(1-cNull(a,1)));
     





%--------------------------------------------------------------------------
% function [vi, d, Ct] = colocalMeasure(imageA, threshold, method)

%Example: Generate Artificial Data to calculate colocalization
% imageA = number of query points
% threshold = distance threshold for colocalization
% method = 1-evenly distributed points, 2- randomly placed points
%          3 - Clustered points

% %Creates imageA number of random points within range, default = 50
%  QP = randi(50,imageA,2);
% % scatter(QP(:,1),QP(:,2),'Marker','.');
% switch method
%     %Evenly distributed query point values-------------------------------------
%     case 1
%         
%         A = [1:50];
%         Y = repmat(A,1,length(A));
%         X = repmat(A,length(A),1);
%         X = X(1:end);
%         X1(:,1) = X;
%         X1(:,2) = Y;
%         %scatter(X1(:,1),X1(:,2),'Marker','.');
%     
% 
%     %--------------------------------------------------------------------------
%     % Random Query Points
%     case 2
%         X1 = randi(50,100,2);
%         %scatter(X1(:,1),X1(:,2),'Marker','.');
%     
%     %--------------------------------------------------------------------------
%     % Clustered Query Point values
%     case 3
%         MU1 = [1 40]; SIGMA1 = [3 0; 0 2];
%         MU2 = [17 21];   SIGMA2 = [3 0; 0 2];
%         MU3 = [25 35];  SIGMA3 = [3 0; 0 3];
%         MU4 = [40 50];  SIGMA4 = [3 0; 0 3];
%         X1 = [mvnrnd(MU1,SIGMA1,20);mvnrnd(MU2,SIGMA2,20);mvnrnd(MU3,SIGMA3,20);
%         mvnrnd(MU4,SIGMA4,20)];
%         scatter(X1(:,1),X1(:,2),'Marker','.');
% end
% 
% 
%     %--------------------------------------------------------------------------
% 
% DT = delaunayTriangulation(X1(:,1),X1(:,2));
% 
% %Find NN
% [vi, d] = nearestNeighbor(DT, QP);
% 
% % Test NN under threshold
% j =1;
% for i = 1:length(d)
%     if d(i) <= threshold
%         C(j) =d(i);
%         j = j+1;
%     end
% end
% 
% %Determine Ct score
% Ct = length(C)/length(d)
% 
% %Produce probablilty density of values
% [f, xi] = ksdensity(d);
% figure()
% plot(xi,f);
% xlim([-5 50])
% ylim([0 1])

 end