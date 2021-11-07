classdef StrainEnergyCalculationProcess < DataProcessingProcess
    % Concrete process for calculating a force field
    %
    % Sangyoon Han & Andrew Jamieson, May 2017
    properties (SetAccess = protected)  
        tMapLimits_
        dELimits_
        distBeadMapLimits_
    end
    
    methods
        function obj = StrainEnergyCalculationProcess(owner,varargin)
            
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = StrainEnergyCalculationProcess.getName;
                super_args{3} = @calculateMovieStrainEnergy;
                if isempty(funParams)
                    funParams=StrainEnergyCalculationProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
        end
        
        function varargout = loadChannelOutput(obj, iFrame, varargin)
            % Input check
            outputList = {'SE_Blobs','SE_FOV','SE_Cell','totalForceBlobs',...
                'totalForceCell','totalForceFOV','mask'};
            ip = inputParser;
            ip.addRequired('obj');
            ip.addRequired('iFrame', @(x) obj.checkFrameNum(x));
            ip.addParameter('useCache',true, @islogical);
            ip.addParameter('output', outputList{3}, @(x) all(ismember(x,outputList)));
            ip.parse(obj, iFrame, varargin{:})
            output = ip.Results.output;
            varargout = cell(numel(output), 1);
            iFrame = ip.Results.iFrame;
            if ischar(output),output={output}; end

            if ~strcmp(output, 'mask')
                % Load all the output data
                s = cached.load(obj.outFilePaths_{1,10}, '-useCache', ip.Results.useCache);

                for iout = 1:numel(output)
                    switch output{iout} 
                        case 'SE_FOV' 
                            val = s.(output{iout}).SE(iFrame);                   
                        case 'SE_Cell'
                            val = s.(output{iout}).SE(iFrame);
                        case 'SE_Blobs'
                            val = s.(output{iout}).SE(iFrame);
                        case 'totalForceBlobs' 
                            val = s.(output{iout}).force(iFrame);                   
                        case 'totalForceCell'
                            val = s.(output{iout})(iFrame);
                        case 'totalForceFOV'
                            val = s.(output{iout})(iFrame);
                        otherwise
                            error('Incorrect Output Var type');
                    end
                    varargout{iout} = val;
                end 
            else
                nFrames = obj.owner_.nFrames_;
                maskFolder = [obj.funParams_.OutputDirectory filesep 'BandMasks'];
                fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
                numStr = @(frame) num2str(frame,fString);
                maskIndPath = @(frame) [maskFolder filesep 'mask' numStr(frame) '.mat'];
                maskObj = load(maskIndPath(iFrame), 'maskOnlyBand');
                varargout{1} = maskObj.maskOnlyBand;
            end
        % @(proc,iChan,iFrame,varargin) proc.getParameters().text;
        % frame = arrayfun(@(x) ['Frame ' num2str(x)], 1:3, 'Uniform' ,0)';
        % t = struct2table(SE_Blobs);
        % t.Properties.RowNames = t.Frame
        % t.Frame = []
        end

        function output = getDrawableOutput(obj)
            i = 1;
            yDist = 20;
            disOpts = {'Color', [1 0 0],...
                       'Position',[10 10],...
                       'FontSize', 14};
            
            output(i).name = 'SE FOV';
            output(i).var = 'SE_FOV';
            output(i).formatData = @(val) struct('String',{{['SE FOV: ' num2str(val)]}}, disOpts{:});
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;
            
            i = i+1;
            disOpts{4} = [10 10+(i-1)*yDist];
            output(i).name = 'SE Cell';
            output(i).var = 'SE_Cell';
            output(i).formatData = @(val) struct('String', {{['SE Cell: ' num2str(val)]}}, disOpts{:});
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;            

            i = i+1;
            disOpts{4} = [10 10+(i-1)*yDist];
            output(i).name = 'SE Blobs';
            output(i).var = 'SE_Blobs';
            output(i).formatData = @(val) struct('String', {{['SE Blobs: ' num2str(val)]}}, disOpts{:});
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;

            i = i+1;
            disOpts{4} = [10 10+(i-1)*yDist];
            output(i).name = 'Total Force FOV';
            output(i).var = 'totalForceFOV';
            output(i).formatData = @(val) struct('String',{{['Total Force FOV: ' num2str(val)]}}, disOpts{:});
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;
            
            i = i+1;
            disOpts{4} = [10 10+(i-1)*yDist];
            output(i).name = 'Total Force Cell';
            output(i).var = 'totalForceCell';
            output(i).formatData = @(val) struct('String', {{['Total Force Cell: ' num2str(val)]}}, disOpts{:});
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;            

            i = i+1;
            disOpts{4} = [10 10+(i-1)*yDist];
            output(i).name = 'Total Force Blobs';
            output(i).var = 'totalForceBlobs';
            output(i).formatData = @(val) struct('String', {{['Total Force Blobs: ' num2str(val)]}}, disOpts{:});
            output(i).type = 'movieOverlay';
            output(i).defaultDisplayMethod = @TextDisplay;            
            
            i = i+1;
            output(i).name = 'Cell Periphery Mask';
            output(i).var='mask';            
            output(i).type='movieOverlay';
            output(i).formatData=@StrainEnergyCalculationProcess.getMaskBoundaries;
            colors = hsv(numel(obj.owner_.channels_));
            output(i).defaultDisplayMethod=@(x) LineDisplay('Color',colors(x,:));                                
        end
    end
%         function status = checkChannelOutput(obj,varargin)
%             
%             status = logical(exist(obj.outFilePaths_{1},'file'));
%             
%         end
    methods (Static)
        function name =getName()
            name = 'Strain Energy Calculation';
        end
        function h = GUI()
            h= @strainEnergyCalculationProcessGUI;
        end
    
        function boundaries = getMaskBoundaries(mask)
            % Format mask boundaries in xy-coordinate system
            b=bwboundaries(mask);
            b2 =cellfun(@(x) vertcat(x,[NaN NaN]),b,'Unif',false);
            boundaries =vertcat(b2{:});
            if ~isempty(boundaries)
                boundaries = boundaries(:,2:-1:1);
            else
                boundaries = [NaN NaN];
            end
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'strainEnergy'];
            funParams.exportCSV=true;
            funParams.performForceBlobAnalysis=true;
            funParams.useFOV=true;
            funParams.useCellMask=true;
            funParams.subcellularTFM=false;
            funParams.bandWidth=2; % in micron

        end
%         function units = getUnits(varargin)
%             units = 'Traction (Pa)';
%         end
    end
end