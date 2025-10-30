function [fractionFAs,fractionFCs,fractionFCFAs,numG2FAs,numG2FCs,numG2FCFAs,numG2] = analyzeFractionFAsfromG2NAs(MD)
%analyzeFractionFAsfromG2NAs(MD) identify all G2 adhesions then see hwo
%many of them actually mature into FAs or FCs.
%   input
%           MD:             MovieData object that has FocalAdhesionPackage
%                           already run.
%   output
%           fraction:       fraction of numG2FAs to all numG2
%           numG2:          the number of all G2 adhesions
%           numG2FA:        the number of G2 adhesion that actually mature
%                           into FAs
% Sangyoon Han, Nov 2020

%% Identify G2 adhesions
% Get the FA package
faPack = MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Get the classes
iClaProc = MD.getProcessIndex('AdhesionClassificationProcess');
classProc = MD.getProcess(iClaProc);
ChannelIndex = classProc.funParams_.ChannelIndex;
idsClassified = load(classProc.outFilePaths_{4,ChannelIndex});
idGroup2 = idsClassified.idGroup2;
% Get tracks with idGroup2
iAdhProc = faPack.getProcess(7);
tracksG2 = iAdhProc.loadChannelOutput(ChannelIndex,'output','tracksNA','wantFullTrack',true,'idSelected',find(idGroup2));

%% Quantify statistics

numG2 = numel(tracksG2);
%% check if there is FAs
% FC: 3, FA: 4
idG2FAs = arrayfun(@(x) any(x.state==4),tracksG2);
numG2FAs = sum(idG2FAs);
idG2FCs = arrayfun(@(x) any(x.state==3),tracksG2);
numG2FCs = sum(idG2FCs);
numG2FCFAs = sum(idG2FCs | idG2FAs);

fractionFAs = numG2FAs/numG2;
fractionFCs = numG2FCs/numG2;
fractionFCFAs = numG2FCFAs/numG2;

end

