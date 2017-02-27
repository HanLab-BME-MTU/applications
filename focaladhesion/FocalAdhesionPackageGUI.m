function varargout = FocalAdhesionPackageGUI(varargin)
% Launch the GUI for the Focal Adhesion Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%

% Andrew R. Jamieson - Feb 2017

if nargin>0 && isa(varargin{1},'MovieList')
    varargout{1} = packageGUI('FocalAdhesionPackage',...
        [varargin{1}.getMovies{:}],varargin{2:end},'ML',varargin{1});
else
    varargout{1} = packageGUI('FocalAdhesionPackage', varargin{:});
end

end