function varargout = plusTipTrackerPackageGUI(varargin)
% Launch the GUI for the u-Track Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Sebastien Besson 5/2011
%

options = {'packageConstr', @PlusTipTrackerPackage};
if nargin>0 && isa(varargin{1},'MovieList')  
    varargout{1} = packageGUI('TrackingPackage',...
        [varargin{1}.getMovies{:}],varargin{2:end},options{:},'ML',varargin{1});
else
    varargout{1} = packageGUI('TrackingPackage',varargin{:},options{:});
end

end