function analyzeAllDatasets(numClassDataSets,doSave)
% let's say, all the data is in one folder, and the files are named:
% analysisConstrSingleCell_wellA4_03.mat

%read in Stack of bead .mat files:
for c=1:numClassDataSets
    [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
           ['Select First .mat file of dataset class',num2str(c)]);

       if ~ischar(filename) || ~ischar(pathname)
           return;
       end

       class(c).inputFileList = getFileStackNames([pathname filesep filename]);
       class(c).numFiles=numel(class(c).inputFileList);
       
       pattern='well';
       startIndex=strfind(filename, pattern);
       lengthID=2; %e.g. A4
       class(c).wellClass=filename((startIndex:startIndex+length(pattern)+lengthID-1));       
end

for c=1:numClassDataSets
    class(c).elEnergyAllWells=[];
    class(c).maxForceAllWells=[];
    class(c).peakForceOfAllWells=[];
    class(c).elEnergyPerSprAreaAllWells=[];
    for i=1:class(c).numFiles
        %read in the data:
        fileStruct=load(class(c).inputFileList{i});
        
        % get maxForce and elEnergy and MaxForceOfAllFrames
        roundBrkts=false;
        try
            roundBrkts=isfield(fileStruct.constrForceField(1), 'vec');
        end
        
        if roundBrkts
            for j=1:length(fileStruct.constrForceField)
                maxForce(j)=max(sqrt(sum(fileStruct.constrForceField(j).vec.^2,2)));
                if isfield(fileStruct.constrForceField(j), 'stats');
                    display('taken the stats-value')
                    elEnergy(j)=fileStruct.constrForceField(j).stats.elEnergy;
                else
                    elEnergy=fileStruct.elEnergy;
                end
            end
            MaxForceOfAllFrames=max(maxForce); 
        else
            for j=1:length(fileStruct.constrForceField)
                maxForce(j)=max(sqrt(sum(fileStruct.constrForceField{j}.vec.^2,2)));
                if isfield(fileStruct.constrForceField{j}, 'stats');
                    display('taken the stats-value')
                    elEnergy(j)=fileStruct.constrForceField{j}.stats.elEnergy;
                else
                    elEnergy=fileStruct.elEnergy;
                end
            end
            MaxForceOfAllFrames=max(maxForce); 
        end          
        
        %collect elastic energy:        
        %elEnergy=fileStruct.elEnergy;
        class(c).elEnergyAllWells=horzcat(class(c).elEnergyAllWells,elEnergy);
        
        %collect maximum force in each frame:
        %maxForce=fileStruct.maxForce;
        class(c).maxForceAllWells=horzcat(class(c).maxForceAllWells,maxForce);
        
        %collect peak force of all frames:
        %MaxForceOfAllFrames=fileStruct.MaxForceOfAllFrames;
        class(c).peakForceOfAllWells=horzcat(class(c).peakForceOfAllWells,MaxForceOfAllFrames);
        
        %collect elastic energy per spreading area for all frames:
        roiExists=false;
        try
            roiExists=isfield(fileStruct.constrForceField{1}, 'ellipse');
        end
        if roiExists
            for j=1:length(fileStruct.constrForceField)
                cellArea(j)=sum(sum(fileStruct.constrForceField{j}.ellipse.mask));
            end
            elEnergyPerSprArea=elEnergy./cellArea;            
            class(c).elEnergyPerSprAreaAllWells=horzcat(class(c).elEnergyPerSprAreaAllWells,elEnergyPerSprArea);
        end
        clear cellArea elEnergy maxForce;
    end
end

%organize the data such that it can be used for a boxplot:
elEnergyBoxPlotVec=[];
maxForceBoxPlotVec=[];
peakForceBoxPlotVec=[];
elEnPerSprArPlotVec=[];

boxPlotGroupVec=[];
boxPlotGroupVecPeakForce=[];

for c=1:numClassDataSets
    %arrange the values for the elastic energy:
    elEnergyBoxPlotVec=horzcat(elEnergyBoxPlotVec,class(c).elEnergyAllWells);
    
    %arrange the values for the maximum force:
    maxForceBoxPlotVec=horzcat(maxForceBoxPlotVec,class(c).maxForceAllWells);
    
    %arrange the values for the peak force:
    peakForceBoxPlotVec=horzcat(peakForceBoxPlotVec,class(c).peakForceOfAllWells);
    
    if roiExists
        %arrange the values for the maximum force:
        elEnPerSprArPlotVec=horzcat(elEnPerSprArPlotVec,class(c).elEnergyPerSprAreaAllWells);
    end 
    
    %generate the grouping vector for the boxplot of el. energy or max force or elastic energy per spreading area:
    boxPlotGroupVec=horzcat(boxPlotGroupVec,repmat({[class(c).wellClass,' (',num2str(class(c).numFiles),')']},1,length(class(c).elEnergyAllWells)));
    
    %generate the grouping vector for the boxplot of the peak force:
    boxPlotGroupVecPeakForce=horzcat(boxPlotGroupVecPeakForce,repmat({[class(c).wellClass,' (',num2str(class(c).numFiles),')']},1,length(class(c).peakForceOfAllWells)));
    
    if nargin>1 && doSave==1
        if ~isdir('csvFiles')
            mkdir('csvFiles')
        end
        % write a csv file for each class:
        % elastic energy:
        % create the filename:
        filename=['csvFiles',filesep,class(c).wellClass,'_elE_','Neq',num2str(class(c).numFiles),'.dat'];
        csvwrite(filename,(class(c).elEnergyAllWells)');

        % max force:
        % create the filename:
        filename=['csvFiles',filesep,class(c).wellClass,'_maxF_','Neq',num2str(class(c).numFiles),'.dat'];
        csvwrite(filename,(class(c).maxForceAllWells)');

        % peak force:
        % create the filename:
        filename=['csvFiles',filesep,class(c).wellClass,'_peakF_','Neq',num2str(class(c).numFiles),'.dat'];
        csvwrite(filename,(class(c).peakForceOfAllWells)');
    else
        display('If you want to save the results in .dat-files use the command: analyzeAllDatasets(n,1)');
    end
end


% for c=1:numClassDataSets
%     figure(1)
%     errorbar(c,mean(class(c).elEnergyAllWells),std(class(c).elEnergyAllWells),'x')
%     hold on
%     title('Elastic energy provided by the cells mean+-std')
% end
% hold off
% 
% for c=1:numClassDataSets
%     figure(2)
%     errorbar(c,mean(class(c).maxForceAllWells),std(class(c).maxForceAllWells),'x')
%     hold on
%     title('Maximum local stress built up by the cells mean+-std')
% end
% hold off
% 
% for c=1:numClassDataSets
%     figure(3)
%     errorbar(c,mean(class(c).peakForceOfAllWells),std(class(c).peakForceOfAllWells),'x')
%     hold on
%     title('Peak stress built up by the cells mean+-std')
% end
% hold off
% 
% if roiExists
%     for c=1:numClassDataSets
%         figure(4)
%         errorbar(c,mean(class(c).elEnergyPerSprAreaAllWells),std(class(c).elEnergyPerSprAreaAllWells),'x')
%         hold on
%         title('Elastic energy per spreading area, mean+-std')
%     end
%     hold off
% end


figure(10)
boxplot(elEnergyBoxPlotVec,boxPlotGroupVec,'notch','marker');
title('Elastic energy provided by the cells')

figure(20)
boxplot(maxForceBoxPlotVec,boxPlotGroupVec,'notch','marker');
title('Maximum local stress built up by the cells')

figure(30)
boxplot(peakForceBoxPlotVec,boxPlotGroupVecPeakForce,'notch','marker');
title('Peak stress built up by the cells')

if roiExists
    figure(40)
    boxplot(elEnPerSprArPlotVec,boxPlotGroupVec,'notch','marker');
    title('Elastic energy per spreading area')
end
