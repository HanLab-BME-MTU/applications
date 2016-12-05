function mapXcorrCurvePermutation_Vel(MD, iChan1, chan1Name, layerMax, figuresDir, varargin) 
% mapXcorrCurvePermutation_Vel Perform cross correlation analysis between a
% channel and the edge velocity. It plots cross correlation maps, their mean curves
% at each layer together with confidence bounds based on permutation, and a
% topograph of the cross correlations at lag 0. The cross correlations at
% lag h are Corr(chan1_{t+h}, Vel_t).
% It computes cross correlations in a fashion that can handle many NaN's 
% by utilizing nanXcorrMaps.m function.
%
% Usage:
%       mapXcorrCurvePermutation_Vel(MD, 1, 'Actin', 4, figuresDir, 'impute', 0, 'parpoolNum', 16) 
%
% Input:
%       MD          - a movieData object
%       iChan1      - a channel index
%       chan1Name   - a short name for the channel. eg. 'Actin'
%       layerMax    - maximum layer to be analyzed      
%       figuresDir  - a directory where plots are saved as png files
%
% Output: png files are saved in the figuresDir.
%       
% Option:
%       figFlag     - if 'on', matlab figures are ploted. Default is 'off'.
%       lagMax      - the maximum lag to compute. Default is the number of
%                   time frames devided by 4.
%       fullRange   - if true, smoothed xcorr maps are given in [-1, 1]
%                   scale to compare multiple xcorr maps. Default is false.
%       impute      - if true, moderate missing values are imputed by using
%                   knnimpute.m function. Default is true.
%       parpoolNum  - number of local parallel pool used during permutation. Default is 4.
%       rseed       - input for running rng('default'); rng(rseed). Default
%                   is 'shuffle'. If it is a specific number, the permutation will give
%                   the same result.
%       numPerm     - number of permutation. Default is 1000.
%       WithN       - if true, it uses an alternative windowSampling result
%                   which is obtained by sampleMovieWindowsWithN.m and includes number
%                   of pixels for each windows. Default is false.
%   
%
% Jungsik Noh, 2016/10/22

tmax = MD.nFrames_;

ip = inputParser; 
ip.addParameter('figFlag', 'off');
ip.addParameter('lagMax', round(tmax/4), @isnumeric);
ip.addParameter('fullRange', false);
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);
ip.addParameter('impute', true);
ip.addParameter('WithN', false);

parse(ip, varargin{:})
p = ip.Results;

%figFlag = p.figFlag;

%%  figuresDir setup
if ~isdir(figuresDir); mkdir(figuresDir); end

%%  getting Maps from channel & vel (ch0)

[~, ~,MDtimeInterval_, wmax, tmax, ~, ~, imActmap1] ...
            = mapOutlierImputation(MD, iChan1, layerMax, 'impute', p.impute, 'WithN', p.WithN); 

[~, ~, ~, ~, ~, ~, ~, imVelmap] ...
            = mapOutlierImputation(MD, 0, 1, 'impute', p.impute); 

%%  adjust actmap according to velmap. Match layers

for indL = 1:layerMax
    imActmap1{indL}(:, 1) = [];             % Start from frame=2
end

imActmap2 = cell(1, layerMax);
for indL = 1:layerMax
    imActmap2{indL} = imVelmap{1}(:, 2:tmax);
end

       

%%  input prepare and call xcorrCurvePermutationTest

ch1Actmap = imActmap1;
ch2Actmap = imActmap2;
ch1ActmapName = chan1Name;
ch2ActmapName = 'Vel';
%lagMax
fsaveName0 = ['Ch', num2str(iChan1), 'Ch', num2str(0)];

xcorrMat = xcorrCurvePermutationTest(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
              fsaveName0, MDtimeInterval_, figuresDir, ...
              'figFlag', p.figFlag, 'fullRange', p.fullRange, 'lagMax', p.lagMax, ...
              'numPerm', p.numPerm, 'parpoolNum', p.parpoolNum, 'rseed', p.rseed);
  


%%  Topographs of xcorr

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);

nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
topoMap = nan(wmax, nBandMax_);

for indL = 1:layerMax
    tmp = xcorrMat{indL};
    tmp2 = tmp(:, p.lagMax+1);                % 7 = lagMax+1+h. lag=gef_t+..., stdN_t
    topoMap(:, indL) = tmp2;
end


title0 = ['xcorr(', ch1ActmapName, '_{t}, ', ch2ActmapName, '_t)'];   % lag = {t+ ...} 
topoFig_xcorrChan1Chan2 = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);
 

%%
saveas(topoFig_xcorrChan1Chan2, fullfile(figuresDir, ['/topoFig_xcorr_', fsaveName0, '.png']), 'png')  


%%
disp('====End of mapXcorrCurvePermutation_Vel====')


end




