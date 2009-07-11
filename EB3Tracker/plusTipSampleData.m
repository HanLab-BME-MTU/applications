function [discrimMatrices]=plusTipSampleData(movGroup,prop2sample,discrimMatrices,doPlot)

switch prop2sample
    case {'growthSpeeds','shrinkSpeeds','lifeTimes','gapLength'} % average KS (mean subtracted) test on pop, t-test on sample means

        nReps=100; % do 100 repetitions
        maxSampleSize=400; % max sample size is 400 for KS test
        nGroups=length(movGroup); % number of groups to sample

        % nGroups-vector of population sizes
        popSize=arrayfun(@(x)length(x.(prop2sample)),movGroup);

        % sample each group with smallest distribution size
        sampSizeLargest=min(popSize);
        sampleDataLargest=arrayfun(@(x)randsample(x.(prop2sample),sampSizeLargest),movGroup,'UniformOutput', false);

        if ~isempty(doPlot)
            % make boxplot of group distributions
            figure
            boxplot(cell2mat(sampleDataLargest'))
            xlabel(' ');
            ylabel(prop2sample);
            set(gca,'XTickLabel',{movGroup.label}');
            rotateticklabel(gca,30);
            title([movGroup(1,1).common2group ' - ' prop2sample]);
            
            plotDir=[doPlot filesep 'boxPlots'];
            if ~isdir(plotDir)
                mkdir(plotDir)
            end
            saveas(gcf,[plotDir filesep movGroup(1,1).common2group '_' prop2sample '.fig']);
            close(gcf)
        end

        % sample each group for maxSampleSize values
        sampSizeKS=min([sampSizeLargest; maxSampleSize]);

        %initialize discrimination matrix
        tempDiscrimMat=zeros(nGroups,nGroups,nReps);

        % KS test with mean subtracted = 11
        testStructure.(prop2sample)=[11 11];

        sMeans=zeros(nReps,nGroups);
        for iRep=1:nReps
            % get maxSampleSize values from each group, put in structure, and find
            % mean to be used in later test
            sampleDataKS=arrayfun(@(x)randsample(x.(prop2sample),sampSizeKS),movGroup,'UniformOutput', false);
            sData=cell2struct(sampleDataKS,prop2sample,2);
            sMeans(iRep,:)=mean([sData.(prop2sample)]);

            % make discrimination matrix for property - KS test of sampled population
            compMatrices=discriminationMatrix(sData,testStructure);
            tempDiscrimMat(:,:,iRep)=compMatrices.(prop2sample);
        end

        M=diag(ones(nGroups,1)); M=swapMaskValues(M,[0 1],[1 nan]);
        lowerTri=tril(ones(nGroups));
        upperTri=triu(ones(nGroups));

        % get average discrimination matrix for property from nRep trials
        discrimMatValues=mean(tempDiscrimMat,3).*M;

        clear testStructure
        clear sData

        % make a structure containing mean values from the nReps trials
        sData=cell2struct(mat2cell(sMeans,nReps,ones(nGroups,1))',[prop2sample 'Means'],2);

        % make discrimination matrix for property - t-test of sampled means
        testStructure.([prop2sample 'Means'])=[1 1];
        tempDiscrimMat=discriminationMatrix(sData,testStructure);
        discrimMatMeans=tempDiscrimMat.([prop2sample 'Means']).*M;

        % store disc matrices in structure
        discrimMatrices.(prop2sample)=discrimMatValues;
        discrimMatrices.([prop2sample 'Means'])=discrimMatMeans;
        
        % combine the two matrices into one, with labels
        tempMerge=num2cell(discrimMatValues.*upperTri + discrimMatMeans.*lowerTri);
        idx=find(eye(size(tempMerge,1)));
        tempMerge(idx)=cellfun(@(x)num2str(x),tempMerge(idx),'UniformOutput',false);
        
        cLabel={movGroup.label};
        rLabel=[{prop2sample}; {movGroup.label}'];
        tempMerge=[rLabel [cLabel; tempMerge]];
        discrimMatrices.([prop2sample 'TKS'])=tempMerge;
        
    otherwise
        error('prop2sample not supported')
end