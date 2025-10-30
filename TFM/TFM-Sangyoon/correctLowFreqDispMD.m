function [] = correctLowFreqDispMD(MD)
% function [] = correctLowFreqDispMD(MD) runs correctLowFreqSigFromField
% and overwrites displacement field

%% Get the displacement field
%% Load the forcefield
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));
%% Load the displfield
iCorrectedDisplFieldProc = 3;
CorrectedDisplFieldProc=TFMPackage.processes_{iCorrectedDisplFieldProc};
if ~isempty(CorrectedDisplFieldProc)
    dispProc = CorrectedDisplFieldProc;
else
    disp('No iCorrectedDisplFieldProc found. Trying to use displacementField from step 2...')
    iCalculatedDisplFieldProc = 2;
    dispProc=TFMPackage.processes_{iCalculatedDisplFieldProc};
end

%% For multiple frame
nFrames = MD.nFrames_;
for ii=1:nFrames
    % Make sure if each tmap has its contents
    curDField=dispProc.loadChannelOutput('output','displField','iFrame',ii);
    curCorrectedField = correctLowFreqSigFromField(curDField);
    
end
% Save backup of the original field

% Overwrite new one
end