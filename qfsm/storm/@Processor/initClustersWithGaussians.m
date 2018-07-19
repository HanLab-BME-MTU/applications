% function initClustersWithGaussians(obj)
% 
% % Loop through all the edges until convergence
% while 1
%     
%     for e=1:size(obj.data.edges,1)
%         
%         e1 = obj.data.edges(e,1); % Cluster index of cluster 1
%         e2 = obj.data.edges(e,2);
%         
%         c1 = obj.data.clusters{e1}; % First subcluster
%         c2 = obj.data.clusters{e2};
%         
%         p1 = obj.data.points(c1,:);
%         p2 = obj.data.points(c2,:);
%         p12 = [p1;p2];
%         
%         w1 = 1./obj.data.error(c1,:);
%         w2 = 1./obj.data.error(c2,:);
%         w12 = [w1;w2];
%         
%         % Fit the Gaussian models
%         maxCurvature = 1; % Dummy curvature
%         [cP1 t1 res1] = TLSFitBezierWeightedConstrainedCP(p1,w1,0,maxCurvature);
%         [cP2 t2 res2] = TLSFitBezierWeightedConstrainedCP(p2,w2,0,maxCurvature);
%         [cP12 t12 res12] = TLSFitBezierWeightedConstrainedCP(p12,w12,0,maxCurvature);
%         
%         % Model complexity (model1,model2,mixture,variance)
%         cpSplit = 3+3+2+2; 
%         cpJoined = 3+0+1;
%         
%         % Compute the likelihood of the models
%         
%                     sumWeightedSqDist = cellfun(@(r) sum(sum((r.*pointWeights).^2),2),res);
%             sumWeightedSqDistC1 = sum(sum((resC1.*pointWeights1).^2),2);
%             sumWeightedSqDistC2 = sum(sum((resC2.*pointWeights2).^2),2);
%             
%             varFun = @(d) (d+2*betaVar)/(2*(betaVar-modeVar)/modeVar+nPoints+2);
%             varEst = arrayfun(varFun,sumWeightedSqDist);
%             varEstC1 = (sumWeightedSqDistC1+2*betaVar)/(2*(betaVar-modeVar)/modeVar+nPoints1+2);
%             varEstC2 = (sumWeightedSqDistC2+2*betaVar)/(2*(betaVar-modeVar)/modeVar+nPoints2+2);
%             
%             fun = @(a,b,c) -nPoints*log((2*pi*c)^0.5*a)-0.5/c*b;
%             logL = arrayfun(fun,length,sumWeightedSqDist,varEst); 
%             logL1 = -nPoints1*log((2*pi*varEstC1)^0.5*lengthC1)-0.5/varEstC1*sumWeightedSqDistC1;
%             logL2 = -nPoints2*log((2*pi*varEstC2)^0.5*lengthC2)-0.5/varEstC2*sumWeightedSqDistC2;
%             logL(maxDegreeBezier+1) = nPoints1*log(nPoints1/nPoints)+logL1+nPoints2*log(nPoints2/nPoints)+logL2;
%             logL = reshape(logL,maxDegreeBezier+1,1);
%         
%         % Compute BIC
%         
%         % Compute the edge weight
%         obj.data.weights = 
%         
%     end
%     
%     % Check for convergence
%     if all(obj.data.weights<0)
%         break;
%     end
%     
%     % Find the maximum matching of the edges
%     matching = maxWeightedMatching(obj.data.nClusters,obj.data.edges,obj.data.weights);
%     matching(obj.data.weights<0) = false;
%     
%     % Merge the clusters
%     obj.data.edges = obj.data.edges(matching,:);
%     obj.data.weights = obj.data.weights(matching,:);
%     
%     % Update the edges
%     
% end
% 
% % For all clusters containing more than 3 points fit a linear model
% 
% 
% end
% 
