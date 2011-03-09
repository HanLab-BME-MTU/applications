classdef BiosensorsPackage < Package
% A concrete process for Biosensor Package
    
methods (Access = public)
    function obj = BiosensorsPackage (owner,outputDir)
           % Construntor of class MaskProcess
           if nargin == 0
              super_args = {};
           else
               % Owner: MovieData object
               super_args{1} = owner;
               super_args{2} = 'Biosensors'; 
               % Dependency Matrix (same length as process class name
               % string)
               super_args{3} = BiosensorsPackage.getDependencyMatrix;
                                
               % Process CLASS NAME string (same length as dependency matrix)
               % Must be accurate process class name
               super_args{4} = {'ThresholdProcess',... 
                                'BackgroundMasksProcess',... 
                                'MaskRefinementProcess',... 
                                'DarkCurrentCorrectionProcess',...
                                'ShadeCorrectionProcess',...
                                'BackgroundSubtractionProcess',...
                                'TransformationProcess',...
                                'BleedthroughCorrectionProcess',...
                                'RatioProcess',...    
                                'PhotobleachCorrectionProcess',...
                                'OutputRatioProcess'...
                                };
                            
               super_args{5} = [outputDir  filesep 'BiosensorsPackage']; 
                
           end
           % Call the superclass constructor 
           obj = obj@Package(super_args{:});
    end
    
    function processExceptions = sanityCheck(obj,full,procID) % throws Exception Cell Array
        % Sanity Check
        % full package panity check: true or false
        nProcesses = length(obj.processClassNames_);
            
        if nargin < 2
            full = true;
            procID = 1:nProcesses;
        end
            
        if nargin < 3
           procID = 1:nProcesses ;
        end
            
        if strcmp(procID,'all')
            procID = 1:nProcesses;
        end
            
        if any(procID > nProcesses)
            error('User-defined: process id exceeds number of processes');
        end
            
        processExceptions = obj.checkProcesses(full,procID);  % throws Exception Cell Array
            
        if full
            
            for i = procID
                if isempty(obj.processes_{i})
                    continue;
                else
                    parentIndex = find(obj.depMatrix_(i,:));
                    
                    if length(parentIndex) > 1
                        
                        switch i
                            case 6
                                parentIndex = 5;
                            case 7
                                parentIndex = 6;
                            case 8
                                parentIndex = 6;
                            case 9
                                parentIndex = 6;
                            case 11
                                parentIndex = 9;
                            otherwise
                                parentIndex = parentIndex(1);
                        end

                    end
                    
                    % Check if input channels are included in dependent
                    % processes
                    
                    if  ~isempty(parentIndex) && ~isempty(obj.processes_{parentIndex})
                        
                        tmp =setdiff(obj.processes_{i}.funParams_.ChannelIndex, ...
                                union(obj.processes_{parentIndex}.funParams_.ChannelIndex,...
                                      find(obj.processes_{parentIndex}.checkChannelOutput)));                            
                                
                        if  ~isempty(tmp)
                            
                            if length(tmp) ==1

                                ME = MException('lccb:input:fatal',...
                                    'Input channel ''%s'' is not included in step %d. Please include this channel in %d step or change another input channel for the current step.',...
                                    obj.owner_.channels_(tmp).channelPath_, parentIndex, parentIndex);
                            else
                                
                                ME = MException('lccb:input:fatal',...
                                    'More than one input channels are not included in step %d. Please include these channels in %d step or change other input channels for the current step.',...
                                    parentIndex, parentIndex);
                            end
%                             processExceptions{i} = horzcat(processExceptions{i}, ME); 
                            processExceptions{i} = [ME, processExceptions{i}];
                            
                        end
                    end
                    
                    % Check the validity of mask channels in Background
                    % Subtraction Process (step 6 in biosensors package)
                    if i == 6
                        if ~isempty(obj.processes_{2})
                            tmp = setdiff(obj.processes_{i}.funParams_.MaskChannelIndex, ...
                                union(obj.processes_{2}.funParams_.ChannelIndex,...
                                      find(obj.processes_{2}.checkChannelOutput)));  
                            
                            if  ~isempty(tmp)
                                if length(tmp) ==1
                                    ME = MException('lccb:input:fatal',...
                                        'The mask channel ''%s'' is not included in step 2 (Background Mask Creation). Please include this channel in step 2 or change to another channel that has mask.',...
                                        obj.owner_.channels_(tmp).channelPath_);
                                else
                                    ME = MException('lccb:input:fatal',...
                                        'More than one mask channels are not included in step 2 (Background Mask Creation). Please include these channels in step 2 or change to other channels that have masks.');
                                end
%                               processExceptions{i} = horzcat(processExceptions{i}, ME); 
                                processExceptions{i} = [ME, processExceptions{i}];
                            end                            
                        end
                    end 
                    
                    % Check the validity of bleed channels in Bleedthrough
                    % Correction Process (step 8 in biosensors package)
                    if i == 8
                        if ~isempty(obj.processes_{6})
                            tmp = setdiff(obj.processes_{i}.funParams_.BleedChannelIndex, ...
                                union(obj.processes_{6}.funParams_.ChannelIndex,...
                                      find(obj.processes_{6}.checkChannelOutput)));  
                            
                            if  ~isempty(tmp)
                                if length(tmp) ==1
                                    ME = MException('lccb:input:fatal',...
                                        'The bleedthrough channel ''%s'' is not included in step 6 (Background Subtraction). Please include this channel in step 6 or change to another bleedthrough channel that is background-subtracted.',...
                                        obj.owner_.channels_(tmp).channelPath_);
                                else
                                    ME = MException('lccb:input:fatal',...
                                        'More than one bleedthrough channels are not included in step 6 (Background Subtraction). Please include these channels in step 6 or change to other bleedthrough channels that are background-subtracted.');
                                end
%                                 processExceptions{i} = horzcat(processExceptions{i}, ME); 
                                processExceptions{i} = [ME, processExceptions{i}];
                            end                            
                        end
                    end 
                    
                    % Check the validity of mask channesl in Ratio Process 
                    % (step 9 in biosensors package)
                    if i == 9
                        if ~isempty(obj.processes_{1})
                            tmp = setdiff(obj.processes_{i}.funParams_.MaskChannelIndex, ...
                                union(obj.processes_{1}.funParams_.ChannelIndex,...
                                      find(obj.processes_{1}.checkChannelOutput)));  
                            
                            if  ~isempty(tmp)
                                if length(tmp) ==1
                                    ME = MException('lccb:input:fatal',...
                                        'The mask channel ''%s'' is not included in step 2 (Background Mask Creation). Please include this channel in step 2 or change to another channel that has mask.',...
                                        obj.owner_.channels_(tmp).channelPath_);
                                else
                                    ME = MException('lccb:input:fatal',...
                                        'More than one mask channels are not included in step 2 (Background Mask Creation). Please include these channels in step 2 or change to other channels that have masks.');
                                end
%                                 processExceptions{i} = horzcat(processExceptions{i}, ME); 
                                processExceptions{i} = [ME, processExceptions{i}];
                            end                            
                        end
                    end
                    
                    % Photobleach and Output step:
                    % Check if input channel (single channel) is the numerator of ratio channel 
                    if i == 10 || i == 11
                        if ~isempty(obj.processes_{9})
                            
                            % ok = [0 1 0 1...] 1: the channel is numerator
                            ok = obj.processes_{9}.checkChannelOutput;
                            
                            if obj.processes_{i}.funParams_.ChannelIndex ~= ...
                                    obj.processes_{9}.funParams_.ChannelIndex(1) && ...
                                    ~ok( obj.processes_{i}.funParams_.ChannelIndex )
                                
                                    ME = MException('lccb:input:fatal',...
                                    'The input channel of current step must be the numerator of ratio channels. There can be multiple numerator channels generated by Ratioing step (step 9) in multiple times of processing.');                                
                                    processExceptions{i} = [ME, processExceptions{i}];
                                
                            end
                        end                        
                    end                   
                    
                end
                
            end
            
            % Hard-coded, when processing processes 2,3,7,9,  add mask
            % process to the processes' funParams_.SegProcessIndex
            %
            % If only segmentation process exists:
            %       funParams.SegProcessIndex = [SegmentationProcessIndex]
            %
            % If segmentation and maskrefinement processes both exist:   
            %       funParams.SegProcessIndex = [MaskrefinementProcessIndex,  SegmentationProcessIndex]
            %

      
            for i = intersect(procID, [2 3 7 9])
                
                if ~isempty(obj.processes_{i})
                    
                    if ~isempty(obj.processes_{1}) % Threshold process
                        
                         funParams = obj.processes_{i}.funParams_;
                         segPI = find(cellfun(@(x)isequal(x, obj.processes_{1}), obj.owner_.processes_));
                         if length(segPI) > 1
                             error('User-defined: More than one identical Threshold processes exists in movie data''s process list.')
                         end
                         funParams.SegProcessIndex = segPI;
                        
                         % If mask transformation or ratioing process, find
                         % if any mask refinement is done
                        if i == 7 || i == 9
                            if ~isempty(obj.processes_{3}) % MaskRefinement process

                                funParams = obj.processes_{i}.funParams_;
                                segPI = find(cellfun(@(x)isequal(x, obj.processes_{3}), obj.owner_.processes_));
                                if length(segPI) > 1
                                    error('User-defined: More than one identical MaskRefinement processes exists in movie data''s process list.')
                                end
                                funParams.SegProcessIndex = cat(2, segPI, funParams.SegProcessIndex);
                             
                            end  
                        end
                        
                        % if ratioing process, find if there is any mask
                        % refinement process
                        if i == 9
                            segPI = find(cellfun(@(x)isa(x, 'MaskTransformationProcess'), obj.owner_.processes_));
                            
                            if ~isempty(segPI)
                                funParams.SegProcessIndex = cat(2, segPI, funParams.SegProcessIndex);
                            end
                        end
                        
                        obj.processes_{i}.setPara(funParams)
                        
                    else
                        funParams = obj.processes_{i}.funParams_;
                        funParams.SegProcessIndex = [];
                        obj.processes_{i}.setPara(funParams)
                    end
                    
                end
            end
            
        end            
    end
    
    function processExceptions = checkOptionalProcess(obj, procRun, procID)
    % The function check if the successfuly processing of optional
    % processes (with id procID) would affect the update status of
    % decendent processes
        nProcesses = length(obj.processClassNames_);
        processExceptions = cell(1,nProcesses);
        processVisited = false(1,nProcesses); 
        
        for i = procID
            if isempty(obj.processes_{i})
                continue
            else
                [processExceptions, processVisited] = ...
                    obj.dfs_optional(i, procRun, processExceptions, processVisited);
            end
        end
    end
    function [processExceptions, processVisited] = ...
                    dfs_optional(obj, i, procRun, processExceptions, processVisited)
                
        childIndex = find(obj.depMatrix_(:,i)');
        
        
        if isempty(childIndex)
            return
        end
        
        for j=childIndex
            
            if processVisited(j)
                continue
            end
            
            if  ~isempty(obj.processes_{j}) && ...
                ~isempty(setdiff(j, procRun)) && ...
                obj.processes_{j}.success_
            
                obj.processes_{j}.setUpdated (false)
                ME = MException('lccb:depe:warn', ...
                        ['The current step is out of date because one of the optional steps changes the input data of current step. '...
                         'Please run again to update your result.']);   
                processExceptions{j} = horzcat(processExceptions{j}, ME);
            end
            
            processVisited(j) = true;
            
            [processExceptions, processVisited] = ...
                    obj.dfs_optional(j, procRun, processExceptions, processVisited);
        end
    end
end
methods (Static)
    
        function m = getDependencyMatrix()
            % Get dependency matrix
            
               %    1 2 3 4 5 6 7 8 9 10 11
               m = [0 0 0 0 0 0 0 0 0 0 0; %1 MasksProcess
                    1 0 0 0 0 0 0 0 0 0 0; %2 BackgroundMasksProcess
                    1 0 0 0 0 0 0 0 0 0 0; %3 MaskRefinementProcess
                    0 0 0 0 0 0 0 0 0 0 0; %4 DarkCurrentCorrectionProcess
                    0 0 0 0 0 0 0 0 0 0 0; %5 ShadeCorrectionProcess
                    0 1 0 0 1 0 0 0 0 0 0; %6 BackgroundSubtractionProcess
                    0 0 0 0 0 1 0 0 0 0 0; %7 TransformationProcess
                    0 0 0 0 0 1 0 0 0 0 0; %8 BleedthroughCorrectionProcess
                    0 0 0 0 0 1 0 0 0 0 0; %9 RatioProcess
                    0 0 0 0 0 0 0 0 1 0 0; %10PhotobleachCorrectionProcess
                    0 0 0 0 0 0 0 0 1 0 0];%11OutputRatioProcess
        end
        
        function id = getOptionalProcessId()
            % Get the optional process id
            id = [3 7 8 10];
        end
end
    
end

