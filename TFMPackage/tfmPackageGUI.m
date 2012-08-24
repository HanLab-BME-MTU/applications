function varargout = tfmPackageGUI(varargin)
% Launch the GUI for the TFM Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%

% Sebastien Besson, Sep 2011
%

if nargin>0 && isa(varargin{1},'MovieList')
    varargout{1} = packageGUI(@TFMPackage,...
        [varargin{1}.getMovies{:}],varargin{2:end},'ML',varargin{1});
else
    varargout{1} = packageGUI(@TFMPackage,varargin{:});
end

end