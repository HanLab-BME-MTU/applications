function [cT, cNull, cCritical,estimatorM, estimatorC] = colocalMeasurePt2Pt(imageAInfo,imageBInfo, threshold)
%COLOCALMEASUREPT2PT measures fraction of coloclization between two punctate channels 
%Input:
% imageInfoA: Detected positions of molecules in channel which will be assigned a nearest neighbor
% 
%
% imageInfoB: Detected positions of molecules in second channel
% 
%
% Threshold: distance threshold for possible interacting objects
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

[cT,cNull,cCritical,estimatorC,estimatorM]=deal(zeros(length(imageAInfo),1));

% Populate query points and DT points with positions taken from images
% xIndex: Information from detection program which reads x position and
% variance of objects in a single frame
%
% yIndex: Information from detection program which reads y position and
% variance of objects in a single frame

for a = 1:length(imageAInfo)
    
    y = imageSet{a,1};%REMEMBER THIS CHANGES
    image  = strcat('/home2/avega/HMECFyn/Fyn/FynT',y);
    I = imread(image);
    
    [mask, maskList] = calcStateDensity(I);
    xIndex = imageAInfo(a).xCoord(:,1);
    %xVar = imageAInfo(a).xCoord(:,2);
    
    yIndex = imageAInfo(a).yCoord(:,1);
    %yVar = imageAInfo(a).yCoord(:,2);
    
    QP = [yIndex xIndex];
    %scatter(xIndex, yIndex);

    iIndex = imageBInfo(a).xCoord(:,1);
    iVar = imageBInfo(a).xCoord(:,2);
    
    jIndex = imageBInfo(a).yCoord(:,1);
    jVar = imageBInfo(a).yCoord(:,2);
    holder = [jIndex iIndex]; 
    %hold on; scatter(iIndex, jIndex,'b');
%     DT = delaunayTriangulation(iIndex, jIndex);
    
 lia1 = ismember(round(QP),maskList,'rows');   
 lia2 = ismember(round(holder),maskList,'rows');
% Multipling lia (binary vector) by defaultM will take replace coord outside
% boundary with zero, last line removes all zeros from vector

 QP(:,1) = lia1.*QP(:,1);
 QP(:,2) = lia1.*QP(:,2);
 QP( ~any(QP,2), : ) = [];
roundedQP = [round(QP(:,1)) round(QP(:,2))]; 
localMask = zeros(size(I,1),size(I,2));
for i = 1:length(roundedQP)
localMask(roundedQP(i,1),roundedQP(i,2)) = 1;
end
 holder(:,1) = lia2.*holder(:,1);
 holder(:,2) = lia2.*holder(:,2);
 holder( ~any(holder,2), : ) = [];
DT = delaunayTriangulation(holder(:,1), holder(:,2));

    %Find Nearest Neighbor between original channels
    [~, d] = nearestNeighbor(DT, QP);

    % Test NN under threshold
    xi = min(d):0.01:max(d);
    [f] = ksdensity(d,xi);
    %Determine density of distances below threshold to find cNull
    test = find(xi<=threshold);
    cT(a,1) = trapz(xi(1:max(test)),f(1:max(test))); 
    
% %     figure; plot(xi,f);

    %Determine Ct score
    %cT(a,1) = length(C)/length(d);

    %Create another channel to measure q(d)
    [~, d1] = nearestNeighbor(DT, sDensity);
%     A = [1:0.25:512];
%     Y = repmat(A,1,length(A));
%     X = repmat(A,length(A),1);
%     X = X(1:end);
%     X1(:,1) = X;
%     X1(:,2) = Y;
%    [~, d1] = nearestNeighbor(DT, X1);
 
    %Produce probablilty density of values
    xi = min(d1):0.01:max(d1);
    
    [f] = ksdensity(d1,xi);
% %     hold on; plot(xi,f);
% % hold off;
    %Determine density of distances below threshold to find cNull
    test = find(xi<=threshold);
    cNull(a,1) = trapz(xi(1:max(test)),f(1:max(test)));

    % Plot pdf
    % figure()
    % plot(xi,f);
    % xlim([-5 50])
    % ylim([0 1])



    %Calculate estimator
    estimatorM(a,1) = log(cT(a,1)/(1-cT(a,1)))-log(cNull(a,1)/(1-cNull(a,1)));


    % Estimate critical parameters based on cNull and interaction size
    cCritical(a,1) = (binoinv(0.95,length(d),cNull(a,1)))/length(d);
    estimatorC(a,1) = log(cCritical(a,1)/(1-cCritical(a,1)))-log(cNull(a,1)/(1-cNull(a,1)));
     
end




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
%         %scatter(X1(:,1),X1(:,2),'Marker','.');
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