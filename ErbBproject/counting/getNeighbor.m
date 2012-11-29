function nb=getNeighbor(site,step,border,varargin)
%GETNEIGHBOR determines lattice point 'step' away from 'site' with under different boundary conditions
%
%   required input arguments:
%         site -> current lattice site
%         step -> where to go
%       border -> total size of lattice
%
%   optional input arguments:
%       mode -> boundary conditions: 'periodic', 'truncated', default 'periodic'
%
%   output argument:
%       nb -> index of neighboring lattice site
%
%   US, 2012/11/20
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('site',@isscalar);
ip.addRequired('step',@isscalar);
ip.addRequired('border',@isscalar);
ip.addParamValue('mode', 'periodic', @(x) any(strcmpi(x, {'periodic', 'truncated'})));

ip.parse(site,step,border,varargin{:});

border=ip.Results.border;

nb=ip.Results.site+ip.Results.step;

if strcmpi(ip.Results.mode, 'periodic')
    if nb < 1
        nb=nb+border;
    elseif nb > border
        nb=nb-border;
    end
elseif strcmpi(ip.Results.mode,'truncated')
    if nb < 1 || nb > border
        nb=-1;
    end
end   

end