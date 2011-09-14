function varargout = tfmPackageGUI(varargin)
% Launch the GUI for the TFM Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%

% Sebastien Besson, Sep 2011
%

varargout{1} = packageGUI(@TFMPackage,varargin{:});

end