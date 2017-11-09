function [ code , differentParameters ] = parameters2code( parameters, varargin)
%parameters2code Converts a parameters struct to code to facilitate
%documentation and batch script writing
%
% parameters is a structure, perhaps from a Process
% defaultParameters (optional) are the default Process parameters which
% will not be included in the code, or as comments

ip = inputParser;
ip.addRequired('parameters',@isstruct);
ip.addOptional('defaultParameters',struct,@isstruct);
ip.addParamValue('name','params',@ischar);
ip.addParamValue('defaultsInComments',true,@islogical);
ip.addParamValue('exclude',{},@iscellstr);
ip.addParamValue('file','',@ischar);
ip.parse(parameters, varargin{:});

differentParameters = structSubstract(parameters, ip.Results.defaultParameters);

for i=1:length(ip.Results.exclude)
    differentParameters = rmfield(differentParameters, ip.Results.exclude{i});
end

code = gencode( differentParameters, ip.Results.name);
code = strjoin( code ,'\n');

if(ip.Results.defaultsInComments && ~isempty(fields(ip.Results.defaultParameters)))
    comments = strcat({'%    '}, gencode(ip.Results.defaultParameters,ip.Results.name) );
    comments = [ {'%    % Default parameters'} comments {''} {''}];
    comments = strjoin( comments ,'\n');
    code = [ comments code];
end

if(~isempty(ip.Results.file))
    fid = fopen(ip.Results.file,'w');
    fwrite(fid,code);
    fclose(fid);
end



end

