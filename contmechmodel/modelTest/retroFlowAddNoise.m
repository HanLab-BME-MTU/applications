function retroFlowAddNoise(mechDir,simuTrueDir,simuNoiseDir, ...
   rNoiseLevel,aNoiseLevel,numRealizations)
%retroFlowAddNoise: This function adds noise to simulated flow field.
%
% SYNOPSIS: 
%    retroFlowAddNoise(mechDir,simuTrueDir,simuNoiseDir, ...
%       rNoiseLevel,aNoiseLevel,numRealizations);
%
% INPUT:
%    mechDir        : Directory for force reconstruction.
%    simuTrueDir    : Directory where the true (no noise) simulated flow field is 
%                     stored.
%    simuNoiseDir   : Directory where the flow with added noise is stored. It can
%                     be a cell array of directories.
%    rNoiseLevel    : Relative noise level in percentage of the flow vector magnitude.
%    aNoiseLevel    : Absolute noise level in percentage of the median flow speed.
%    numRealizations: Number of realizations of the random noise with the 
%                     specified noise level.
%
%    Note: The length of 'simuNoiseDir', 'rNoiseLevel' and 'aNoiseLevel' has
%          be the same.
%
% OUTPUT:
%    Flow with added noise is saved in 'simuNoiseDir'.
%
% AUTHOR: Lin Ji
% DATE  : Nov. 7, 2006.

if ~iscell(simuNoiseDir)
   simuNoiseDir = {simuNoiseDir};
end

if length(simuNoiseDir) ~= length(rNoiseLevel) || ...
   length(simuNoiseDir) ~= length(aNoiseLevel)
   error(['The length of ''simuNoiseDir'', ' ...
      '''rNoiseLevel'' and ''aNoiseLevel'' has to be the same']);
end

%Load the true simulated flow field.
s = load([mechDir filesep simuTrueDir filesep 'rawDispField' filesep 'rawDispField.mat']);
rawDispField = s.rawDispField;

numNoiseLevels = length(simuNoiseDir);

simIndexForm  = sprintf('%%.%dd',max(length(num2str(numRealizations)),2));

simU   = rawDispField.v;
simSpd = sqrt(sum(simU.^2,2));
medSpd = rawDispField.medSpd;

%Boundary displacement field without noise.
edgD = rawDispField.edgD;

%Add random noise to the simulated flow field.
for kk = 1:numNoiseLevels
   if ~isdir([simuNoiseDir{kk} filesep 'rawDispField'])
      success = mkdir(simuNoiseDir{kk},'rawDispField');
      if ~success
         error('Trouble making directory.')
      end
   end
   sDispFieldDir  = [mechDir filesep simuNoiseDir{kk} filesep 'rawDispField'];

   for ii = 1:numRealizations
      %First, add noise to domain display field.
      % Absolute noise to first component.
      absNoise = medSpd*aNoiseLevel(kk)*randn(size(simU(:,1)));

      % Relative noise to first component.
      relNoise = rNoiseLevel(kk)*simSpd.*randn(size(simU(:,1)));

      rawDispField.v(:,1) = simU(:,1) + absNoise + relNoise;

      %Noise to second component.
      absNoise = medSpd*aNoiseLevel(kk)*randn(size(simU(:,2)));
      relNoise = rNoiseLevel(kk)*simSpd.*randn(size(simU(:,2)));
      rawDispField.v(:,2) = simU(:,2) + absNoise + relNoise;

      %Then, add noise to boundary displacement.
      for jj = 1:length(edgD)
         U1 = edgD(jj).U1;
         U2 = edgD(jj).U2;
         ULen = sqrt(U1.^2+U2.^2);

         absNoise = medSpd*aNoiseLevel(kk)*randn(size(U1));
         relNoise = rNoiseLevel(kk)*ULen.*randn(size(U1));
         U1 = U1 + absNoise + relNoise;

         absNoise = medSpd*aNoiseLevel(kk)*randn(size(U2));
         relNoise = rNoiseLevel(kk)*ULen.*randn(size(U2));
         U2 = U2 + absNoise + relNoise;

         edgD(jj).U1 = U1;
         edgD(jj).U2 = U2;
      end
      rawDispField.edgD = edgD;
      rawDispField.absNoiseLevel = aNoiseLevel(kk);
      rawDispField.relNoiseLevel = rNoiseLevel(kk);

      sDispFieldFile = [sDispFieldDir filesep 'rawDispField' ...
         sprintf(simIndexForm,ii) '.mat'];
      save(sDispFieldFile,'rawDispField');
   end
end


