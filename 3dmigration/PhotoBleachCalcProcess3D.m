classdef PhotoBleachCalcProcess3D < ImageAnalysisProcess
%PHOTOBLEACHCALCPROCESS3D calculates photobleach correction curve for a 3D
%movie, using only the intensities in the masked areas.
%
%
% Hunter Elliott
% 1/2013
%


    properties (SetAccess = protected, GetAccess = public)
        
        
    end
    
    methods (Access = public)
        
        function obj = PhotoBleachCalcProcess3D(owner,funParams)
                                              
            if nargin == 0
                super_args = {};
            else                                                
                
                super_args{1} = owner;
                super_args{2} = 'PhotoBleachCalc';
                super_args{3} = @calc3DMoviePhotobleaching;                               
                
                if nargin < 2 || isempty(funParams)                                       
                    
                    funParams = PhotoBleachCalcProcess3D.getDefaultParams(owner);                 
                                    
                end
                
                super_args{4} = funParams;    
                                
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
        end
        function OK = checkChannelOutput(obj,iChan)
            
           %Checks if the selected channel has valid output files           
           if nargin < 2 || isempty(iChan) || numel(iChan) > 1 || ... 
                   ~all(obj.checkChanNum(iChan));
               error('You must specify a valid channel number!')
           end
           OK =  arrayfun(@(x)(exist(obj.outFilePaths_{x},'file')),iChan);
        end
        
        
        function mg = loadChannelOutput(obj,iChan)
            
            if ~obj.checkChanNum(iChan) || numel(iChan) > 1
                error('You must specify a single, valid channel number!')
            elseif ~obj.checkChannelOutput(iChan)
                error('Specified channel does not have valid photbleaching calc files!')
            end                        
            
            mg = load(obj.outFilePaths_{iChan});
            
            fNames = fieldnames(mg);
            if numel(fNames) ~=1
                error('Invalid photo bleach calculation file!')
            end
            mg = mg.(fNames{1});
        end                
           
    end
    methods(Static)
        function getName
        end        
        function h = GUI()   
            h = [];
        end
        function procClasses = getConcreteClasses()
            procClasses = [];
        end        
        function funParams = getDefaultParams(owner)
              %----Defaults----%      

                    funParams.ChannelIndex = 1:numel(owner.channels_);                                        
                    funParams.OutputDirectory = ...
                        [owner.outputDirectory_  filesep 'photobleach_calc'];                    
                    funParams.BatchMode = false;                                                   
        end
    end            
    
end