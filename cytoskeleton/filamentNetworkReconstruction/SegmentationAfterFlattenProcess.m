classdef SegmentationAfterFlattenProcess < MaskProcess
    % An abstract superclass of all segmentation processes

    % Sebastien Besson 4/2011
    
    methods (Access = protected)
        function obj = SegmentationProcess(owner,name,funName, funParams,...
                outFilePaths)
            % Constructor of class SegmentationAfterFlattenProcess
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            if nargin > 2
                super_args{3} = funName;
            end
            if nargin > 3
                super_args{4} = funParams;
            end
            if nargin > 5
                super_args{5} = outFilePaths;
            end
            % Call the superclass constructor - these values are private
            obj = obj@MaskProcess(super_args{:});
           
        end
    end
    methods(Static)
        function name =getName()
            name = 'Segmentation';
        end
        function h = GUI()
            h= @abstractProcessGUI;
        end
        function procClasses = getConcreteClasses()
            procClasses = ...
                {'ThresholdProcess';
                'MSSSegmentationProcess'};
        end
        
    end
end
