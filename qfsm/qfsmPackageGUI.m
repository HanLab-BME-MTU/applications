function varargout = qfsmPackageGUI(varargin)
% Launch the GUI for the QFSM Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Sebastien Besson 5/2011
%

if nargin>0 && isa(varargin{1},'MovieList')
    varargout{1} = packageGUI('QFSMPackage',...
        [varargin{1}.getMovies{:}],varargin{2:end},'ML',varargin{1});
else
    varargout{1} = packageGUI('QFSMPackage', varargin{:});
end

end