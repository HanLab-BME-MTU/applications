function varargout = FilamentAnalysisPackageGUI(varargin)
% Launch the GUI for the Filament Analysis Package
%
% This function calls the generic packageGUI function, passes all its input
% arguments and returns all output arguments of packageGUI
%
%
% Liya Ding  06/2012

varargout{1} = packageGUI(@FilamentAnalysisPackage,varargin{:});

end