classdef BioSensorsPackage < Package
% A concrete process for BioSensor Package
    
    methods (Access = public)
        function obj = BioSensorsPackage (owner)
           % Construntor of class MaskProcess
           if nargin == 0
              super_args = {};
           else
               % Owner: MovieData object
               super_args{1} = owner;
               super_args{2} = 'BioSensors'; 
               % Dependency Matrix (same length as process class name string)
               super_args{3} = [0 0 0 0
                                1 0 0 0
                                0 1 0 0
                                0 0 1 0];
               % Process CLASS NAME string (same length as dependency matrix)
               % Must have accurate process class name
               super_args{4} = {'SegmentationProcess', ... 
                                'SecondProcess', ...
                                'ThirdProcess', ...
                                'FourthProcess'};
           end
           % Call the supercalss constructor with empty cell array (no
           % argument) if nargin == 0
           obj = obj@Package(super_args{:});
        end
        function procEx = sanityCheck(obj,full,procID) % throws Exception Cell Array
            % Sanity Check
            % full: true or false
            if nargin < 2
                full = true;
            end
            procEx = obj.checkProcesses(full,procID);  % throws Exception Cell Array
        end
    end
    methods (Static)
        function text = getHelp(obj)
            % Function return a string of help text
            text = 'This is specific help text of biosensors package.';
        end
    end
end

