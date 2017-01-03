function [ output_args ] = helperCompareWindTypes( projList, outDir)
% will make a boxplot for all window types run  
% for each movie for a given veil parameter 


 params{1} = 'mednVeloc';
    params{2} = 'persTime';
    
    analType{1} = 'protrusionAnalysis';
    analType{2} = 'retractionAnalysis';


    for iAnal = 1:2
        for iParam =  1:2
             count = 1;
            for iProj= 1:numel(projList)
                
               
                load([projList{iProj} filesep 'GrowthConeAnalyzer' filesep 'movieData.mat']);
                
                
                
                
                % find all windows
                idxWind = find(cellfun(@(x) sum(strcmpi(x.name_,'Windowing')),MD.processes_));
                idxProt = find(cellfun(@(x) sum(strcmpi(x.name_,'Protrusion Sampling')),MD.processes_));
                % make sure for every window there is a corresponding protrusion
                if length(idxWind) == length(idxProt)
                    
                    % search for protrusion samples files
                    for iType = 2:length(idxProt)
                        % load the edgeVelocity quantification mat
                        protDir =  MD.processes_{idxProt(iType)}.funParams_.OutputDirectory;
                        edgeDir = [protDir filesep 'EdgeVelocityQuantification'];
                        load([edgeDir filesep 'EdgeMotion.mat']);
                        dataCellAllProjs{count} = analysisResults.(analType{iAnal}).total.(params{iParam});
                      
                        windProc=  MD.processes_{idxProt(iType)-1};
                        windSize = windProc.funParams_.ParaSize; 
                       
                        windType = windProc.funParams_.MethodName; 
                        reInit = windProc.funParams_.ReInit; 
                        
                        forLabels{count} = ['WindSize ' num2str(windSize) 'Method ' windType 'ReInit' num2str(reInit)]; 
                          count = 1+count; % put them all on the same plot for now
                    end
                    grouping{iProj} = repmat(iProj,length(idxProt)-1,1); % get the grouping idx for each project (each project will be a different color) 
                    
                    
                    
                    
                    % save the edge velocity quantification in the protrusion folder
                    % find all the protrusion sampling processes
                    % find the corresponding windows
                    
                    % load the edge velocity quantification from the MD folder
                    
                    
                else
                    display(['Check ' projList{iProj} 'Number of Windows Run Does NOT match the number of protrusion samples'])
                end
                
                
                
            end
            setAxis('on');
            % just make a large data mat for the boxplot
            forBox = reformatDataCell(dataCellAllProjs);
            groupingAll = vertcat(grouping{:});
            N = size(projList,1); 
            colors = linspecer(N,'sequential'); 
            h1 =  boxplot(forBox,'colorGroup',groupingAll,'notch','on',...
                'outlierSize',1,'colors',colors,'Labels',forLabels,'labelorientation','inline','symbol','+');
            % find the max and min 
            test = arrayfun(@(i) prctile(forBox(:,i),95),1:length(forBox(1,:)));
            xLim = length(forBox(1,:)) +0.5; 
            axis([0.5,xLim,0,max(test)]); 
            set(h1(:),'Linewidth',1.5);
            ylabel({analType{iAnal} ; params{iParam}}); 
            saveas(gcf,[outDir filesep analType{iAnal} '_' params{iParam}]); 
            saveas(gcf,[outDir filesep analType{iAnal} '_' params{iParam} '.tif']); 
            
            
        end
    end



end

