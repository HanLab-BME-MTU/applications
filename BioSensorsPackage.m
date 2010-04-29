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
               super_args{2} = 'BioSensorsPackage'; 
               % Dependency Matrix (same length as process class name string)
               super_args{3} = [0 0 0 0
                                1 0 0 0
                                0 1 0 0
                                0 0 1 0];
               % Process class name string
               super_args{4} = {'Mask', ... 
                                'Protrusion', ...
                                'FinalProcess'};
           end
           % Call the supercalss constructor with empty cell array (no
           % argument) if nargin == 0
           obj = obj@Package(super_args{:});
        end
        function sanityCheck(obj,full) % throws Exception Cell Array
            % Sanity Check
            % full: true or false
            if nargin < 2
                full = false;
            end
            obj.checkProcesses(full)  % throws Exception Cell Array
        end
    end
end

