function varargout = plusTipTrackerPackageGUI(varargin)
% Launch the GUI for the u-Track Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Sebastien Besson 5/2011
%

if isa(varargin{1},'MovieList')
    varargout{1} = packageGUI(@PlusTipTrackerPackage,[varargin{1}.getMovies{:}],varargin{2:end});
else
    varargout{1} = packageGUI(@PlusTipTrackerPackage,varargin{:});
end

end