function varargout = windowingPackageGUI(varargin)
% Launch the GUI for the Windowinf Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%

% Sebastien Besson, July 2011

varargout{1} = packageGUI(@WindowingPackage,varargin{:});

end