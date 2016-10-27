function [corticalData,projData] = mitoticAddEndOnInfoWithClassifications( corticalData,projData,varargin)
% mitoticaddEndOnInfoWithClassifications: adds endOn dwell information to
% the corticalData structure. 
%
% INPUT:
%  % corticalData: (REQUIRED) : structure containing the cortical data
%                               information including the lateral
%                               classification information
% 
%    projData :  (REQUIRED) : output of plusTipSubRoiExtractTracksMITOTICPACKAGE_multMasks.m
%                             contains the information corresponding to the
%                             dwell lifetime and the classification 
%                             (terminal,shrinkage,or pause) of each
%                             subTrack in the region
%    
% % Parameters for classification cut-offs : 
% Note:  if mitoticAddPercentTransitionField has not been previously run for
%        the below input combination of values it will error as it will not be 
%        able to find the lateral classification indices. 
% 
% minDispVect: (PARAM) : numeric
%                        Default : 3 pixels
%                        the min MT Track Length of interest 
%                        (currently in pixels)
%                        see Supp Figure 1 A of Kwon et al. Dev Cell 2015
%
% orientCutOff: (PARAM) : numeric 
%                         Default : 60 (in degrees) 
%                         Local Direction MT Growth relative to the local 
%                         classifications used cut-offs as described in supp 
%                         see Figure S7A of Kwon et al. Dev Cell 2015
%                         
%
% dispInRegion: (PARAM) : numeric 
%                         Default : 0.7 ( in um) 
%                         MT Cort Dist < which track will be automatically be classified
%                         as end-on. 
%                         see Figure S7A of Kwon et al. Dev Cell 2015
%
% dispDiscard : (PARAM) : numeric 
%                         Default : 0.3 (im um)
%                         MT Cort Dist < which the track will be discarded.                         
%                         see Figure S7A of Kwon et al. Dev Cell 2015
% 
% OUTPUT : 
% Adds to corticalData.mat 
%                                           .params(minDispVect_orientCutOff_dispInRegion_dispDiscard)
%                                               .trackInfoEndOnST
%                                                             .dwells 
%                                                                   a rx1 vector where r is
%                                                                   the number of MT subTracks
%                                                                   extracted
%                                                                   in the
%                                                                   given
%                                                                   subRegion
%                                                                   that are classified as
%                                                                   end-on given the
%                                                                   classification conditions
%                                                                   and terminate in either 
%                                                                   a MT terminal or an MT shrinkage
%                                                                   event
%                                               .trackInfoEndOnP
%                                                             .dwells  
%                                                                   a rx1 vector where r is
%                                                                   the number of MT subtracks 
%                                                                   extracted in the given subRegion 
%                                                                   that are classified as
%                                                                   end-on given the
%                                                                   classification conditions
%                                                                   and
%                                                                   terminate
%                                                                   in a
%                                                                   pause.
%% Check INPUT
ip = inputParser;

ip.CaseSensitive = false;

% Parameters for the classification cut-offs
ip.addParameter('minDispVect',3,@(x) isnumeric(x));
ip.addParameter('orientCutOff',60,@(x) isnumeric(x));
ip.addParameter('dispInRegion',0.7,@(x) isnumeric(x));
ip.addParameter('dispDiscard', 0.3,@(x) isnumeric(x));
ip.addParameter('orientDiscard',150, @(x) isnumeric(x)); 

ip.parse(varargin{:});

%% SetUp
minDispVectStr = num2str(ip.Results.minDispVect);
orientStr = num2str(ip.Results.orientCutOff);
dispInRegionStr = num2str(ip.Results.dispInRegion);
dispDiscardStr= num2str(ip.Results.dispDiscard);
myStrdiscardOrient = num2str(ip.Results.orientDiscard);

% change any decimals to 'pt'
dispInRegionStr = strrep(dispInRegionStr,'.','pt');
dispDiscardStr = strrep(dispDiscardStr,'.','pt');
fieldnameC = ['params' minDispVectStr '_' orientStr '_' dispInRegionStr '_' dispDiscardStr '_' myStrdiscardOrient ];
fieldname2C = [minDispVectStr '_' orientStr '_' dispInRegionStr '_' dispDiscardStr '_' myStrdiscardOrient]; % for saving in projData for now
idxLatLog = corticalData.(fieldnameC).idxLatLogical;

dataClass = corticalData.(['OrientVsDispVsDispMTVect_' minDispVectStr]); 
orient = dataClass(:,1); 

% load dwell
dwells = projData.dwellAllTracks;
disp = vertcat(projData.dispByFrame{:});

% extract pause,shrinkage,terminal,undefined...
dataMat = projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix; % load to find which tracks end in term,pause, or shrinkage

dataMat = dataMat(dataMat(:,9) ~= 0,:); % filter the dataMat to get rid of pause, shrinkage, and undefined gaps.

% get the indices for each class (indices are for dataMat filtered)
idxEndTerm = dataMat(:,9) == 1; % collect the terminal growth events
idxEndPause = dataMat(:,9) == 2;
idxEndShrink = dataMat(:,9) == 3; % collect
idxEndUndefined = dataMat(:,9) ==4;


% filter by displacement
dwells = dwells(disp>ip.Results.dispDiscard & (orient<ip.Results.orientDiscard));



if ~isempty(dwells)
    % ok I know this is a bit weird here but filter by these when I do the
    % classifications
    idxEndTerm = idxEndTerm(disp>ip.Results.dispDiscard & orient < ip.Results.orientDiscard);
    idxEndPause = idxEndPause(disp>ip.Results.dispDiscard   & orient < ip.Results.orientDiscard);
    idxEndShrink = idxEndShrink(disp>ip.Results.dispDiscard  & orient < ip.Results.orientDiscard);
    idxEndUndefined = idxEndUndefined(disp>ip.Results.dispDiscard & orient < ip.Results.orientDiscard);
    % get only end on dwells
    dwellsEndOnP = dwells(~idxLatLog & idxEndPause);
    dwellsEndOnST = dwells(~idxLatLog & (idxEndTerm | idxEndShrink));
    nDwellsEndOnP = sum(~idxLatLog & idxEndPause);
    nDwellsEndOnU = sum(~idxLatLog & idxEndUndefined);
else
    dwellsEndOnP = [];
    nDwellsEndOnP = [];
    dwellsEndOnST = [];
    nDwellsEndOnU = [];
end

% record mean params
corticalData.(fieldnameC).stats.EndOn_DwellTimeST_mean = mean(dwellsEndOnST);
corticalData.(fieldnameC).stats.EndOn_DwellTimeP_mean = mean(dwellsEndOnP);
corticalData.(fieldnameC).stats.N_EndOn_DwellTimeP = nDwellsEndOnP;

% for now also save it to the stats
projData.stats.(['EndOn_DwellTimeST' (fieldname2C)]) = mean(dwellsEndOnST);
projData.stats.(['EndOn_DwellTimeP' (fieldname2C)]) = mean(dwellsEndOnP);
%projData.stats.(['PercentPauseAfterEndOn' (fieldname2C)]) = percentPause;
projData.stats.(['N_EndOn_DwellTimeP' (fieldname2C)]) = nDwellsEndOnP;
projData.stats.(['N_EndOn_DwellTimeU' (fieldname2C)]) = nDwellsEndOnU;
% record track params.
corticalData.(fieldnameC).trackInfoEndOnST.dwells = dwellsEndOnST;
corticalData.(fieldnameC).trackInfoEndOnP.dwells = dwellsEndOnP;


end

