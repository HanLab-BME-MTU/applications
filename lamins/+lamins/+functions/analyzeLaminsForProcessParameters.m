function [ ip ] = analyzeLaminsForProcessParameters( ip )
%analyzeLaminsForProcessParamters Add input parser parameters for
%analyzeLaminsForProcess

    if(nargin < 1)
        ip = inputParser;
    end

        ip.addParameter('steerable_sigma',5,@(x) isnumeric(x) && isscalar(x));
        ip.addParameter('channels',[]);
        ip.addParameter('analysisDate','2015_11_25');
        ip.addParameter('outFilePaths',[]);
        ip.addParameter('tz',[]);
        ip.addParameter('clearOutputDir',false,@islogical);
        ip.addParameter('plot',true,@islogical);
        ip.addParameter('output',[]);

end

