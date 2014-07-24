function [meanDist,fragmentation,crosslinking,falsePositivesFraction] = reconstructionScore(obj,sampleInterval,edgeRadius)

disp('Analysis: Computing reconstruction score ...');

tStart = tic;

% Check if the data contains a reference model
if ~isempty(obj.data.simModelBezCP)

    % Check if the data contains a reconstructed model
    if ~isempty(obj.data.modelBezCP)
                
        disp('Analysis: 1/5');
        
        % Sample the reference model
        nRefMod = numel(obj.data.simModelBezCP);
        
        lenRefMod = cellfun(@lengthBezier,obj.data.simModelBezCP);
        nSamples = ceil(lenRefMod/sampleInterval);
        samplePos = arrayfun(@(a) linspace(1/(2*a),1-1/(2*a),a)',nSamples,'UniformOutput',false);
        samples = cellfun(@(a,b) renderBezier(a,arcLengthToNativeBezierParametrization(a,b)),obj.data.simModelBezCP,samplePos,'UniformOutput',false);
        samplesVector = vertcat(samples{:});
        lastSampleIdx = cumsum(nSamples);
        firstSampleIdx = lastSampleIdx-nSamples+1;
        
        disp('Analysis: 2/5');
        
        % Find the neighbor points of the sample points
        radii = repmat(edgeRadius,sum(nSamples),1);
        neighborPntIdx = KDTreeBallQuery(obj.data.points,samplesVector,radii);
        
        disp('Analysis: 3/5');

        % Remove the neighbors which are unclustered points
        neighborPntIdx = cellfun(@(a) setdiff(a,obj.data.nullCluster),neighborPntIdx,'UniformOutput',false);

        % Find the parent model of the neighbor points
        parents = obj.data.parents;
        neighborModIdx = cellfun(@(a) parents(a),neighborPntIdx,'UniformOutput',false);

        % Remove duplicate entries
        neighborModIdx = cellfun(@unique,neighborModIdx,'UniformOutput',false);
        nNeighborModIdx = cellfun(@numel,neighborModIdx);
        neighborModIdxVector = vertcat(neighborModIdx{:});
        lastNeighborModIdx = cumsum(nNeighborModIdx);
        firstNeighborModIdx = lastNeighborModIdx-nNeighborModIdx+1;

        % Sample indices corresponding to neighborModIdxVector
        sampleIdx = arrayfun(@(a,b) ones(b,1)*a,(1:sum(nSamples))',nNeighborModIdx,'UniformOutput',false);
        sampleIdxVector = vertcat(sampleIdx{:});
        
        % Compute the false negatives
        % isFalseNegative = cellfun(@isempty,neighborModIdx);

        disp('Analysis: 4/5');
        
        % Compute the distance between the sample points and all the
        % neighbor models           
        [distVector,tVector] = arrayfun(@(a,b) distancePointBezier(obj.data.modelBezCP{b},samplesVector(a,:)),sampleIdxVector,neighborModIdxVector);
                
        % Reconstruct dist cell array
        dist = arrayfun(@(a,b) distVector(a:b),firstNeighborModIdx,lastNeighborModIdx,'UniformOutput',false);
        
        % Find the closest model
        [distMin,closestMod] = cellfun(@min,dist,'UniformOutput',false);
        isNotWithinReach = cellfun(@isempty,distMin);
        
        closestModIdxVector = cell(numel(closestMod),1);
        closestModIdxVector(~isNotWithinReach) = cellfun(@(a,b) a(b),neighborModIdx(~isNotWithinReach),closestMod(~isNotWithinReach),'UniformOutput',false);
            
        % Reconstruct closestModIdx cell array
        closestModIdx = arrayfun(@(a,b) closestModIdxVector(a:b),firstSampleIdx,lastSampleIdx,'UniformOutput',false);
        
        % Warn if sample points do not have neighbor models closer than edgeRadius
        nNotWithinReach = nnz(isNotWithinReach);
        if nNotWithinReach 
            fprintf('Analysis: %d/%d samples have no model within reach!\n',nNotWithinReach,sum(nSamples));
        end
        
        % Remove empty fields of the inner cell array
        closestModIdx = cellfun(@(a) a(~cellfun(@isempty,a)),closestModIdx,'UniformOutput',false);
        
        % Remove empty fields of the outer cell array
        emptyRefModFields = cellfun(@isempty,closestModIdx);
        closestModIdx = closestModIdx(~emptyRefModFields);
        refModIdx = (1:nRefMod)'; refModIdx = refModIdx(~emptyRefModFields);
        
        % Convert inner cell array to an array
        closestModIdx = cellfun(@(a) vertcat(a{:}),closestModIdx,'UniformOutput',false);
        
        disp('Analysis: 5/5');
        
        % Compute the fragmentation index
        [~,~,n] = cellfun(@unique,closestModIdx,'UniformOutput',false);
        counts = cellfun(@(a) accumarray(a(:),1),n,'UniformOutput',false);
        nModels = cellfun(@(a) 1/sum((a/sum(a)).^2),counts);
        fragmentation = mean(nModels);
        disp('=========== RECONSTRUCTION SCORE ===========');
        fprintf('Fragmentation index: %.2f\n',fragmentation);
        
        % Compute crosslinking index
        closestModIdxVector = vertcat(closestModIdx{:});
        nClosestMod = cellfun(@numel,closestModIdx);
        
        groundTruthModelIdx = arrayfun(@(a,b) ones(b,1)*a,refModIdx,nClosestMod,'UniformOutput',false);
        groundTruthModelIdxVector = vertcat(groundTruthModelIdx{:});
        
        nRecMod = numel(obj.data.modelBezCP); 
        rec2GT = cell(nRecMod,1);
        for i=1:numel(closestModIdxVector)
            rec2GT(closestModIdxVector(i)) = {[rec2GT{closestModIdxVector(i)} groundTruthModelIdxVector(i)]};
        end
        
        isFalsePositive = cellfun(@isempty,rec2GT);
        [~,~,n] = cellfun(@unique,rec2GT(~isFalsePositive),'UniformOutput',false);
        counts = cellfun(@(a) accumarray(a(:),1),n,'UniformOutput',false);
        nGT = cellfun(@(a) 1/sum((a/sum(a)).^2),counts);
        crosslinking = mean(nGT);
        fprintf(' Crosslinking index: %.2f\n',crosslinking);

        % Find the false positives (Never the closest model)
        nFalsePositives = nnz(isFalsePositive);
        falsePositivesFraction = nFalsePositives/nRecMod;
        fprintf('    False positives: %.1f %%\n',falsePositivesFraction*100);
        
        % Compute the mean distance for each fragment
        meanDist = mean(cell2mat(distMin(~isNotWithinReach)));
        fprintf('      Mean distance: %.1f\n',meanDist);
        disp('============================================');
        
        %         % Display (Debug)
        %         t = arrayfun(@(a,b) tVector(a:b),firstNeighborModIdx,lastNeighborModIdx,'UniformOutput',false);
        %         tMin = cellfun(@(a,b) a(b),t,num2cell(closestMod));
        %         p = arrayfun(@(a,b) renderBezier(obj.data.modelBezCP{b},a),tMin,closestModIdxVector,'UniformOutput',false);
        %         p = vertcat(p{:});
        %         dis = Show(obj.data);
        %         dis.models();
        %         dis.modelGroundTruth();
        %         dis.imaris.displayPoints(samplesVector,1);
        %         dis.imaris.displayPoints(p,1,[0.0,1.0,1.0,0.0]);
        %         dis.imaris.displaySegments(samplesVector,p,'link');
        
    else
        disp('Analysis: No reconstructed model found!')
    end
else
    disp('Analysis: No reference model found!')
end

tEnd = round(toc(tStart));
fprintf('Analysis: Reconstruction score: Done! %s\n',secs2hms(tEnd));

end






