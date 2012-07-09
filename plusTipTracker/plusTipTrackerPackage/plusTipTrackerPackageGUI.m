function varargout = plusTipTrackerPackageGUI(varargin)
% Launch the GUI for the u-Track Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Sebastien Besson 5/2011
%

options = {'packageName', 'TrackingPackage'};
if nargin>0 && isa(varargin{1},'MovieList')  
    varargout{1} = packageGUI(@PlusTipTrackerPackage,...
        [varargin{1}.getMovies{:}],varargin{2:end},options{:});
else
    varargout{1} = packageGUI(@PlusTipTrackerPackage,varargin{:},options{:});
end

end