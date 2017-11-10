function varargout = windowingPackageGUI(varargin)
% Launch the GUI for the Windowinf Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%

% Sebastien Besson, July 2011

if nargin>0 && isa(varargin{1},'MovieList')
    varargout{1} = packageGUI('WindowingPackage',...
        [varargin{1}.getMovies{:}],'ML',varargin{1},varargin{2:end});
else
    varargout{1} = packageGUI('WindowingPackage',varargin{:});
end
end