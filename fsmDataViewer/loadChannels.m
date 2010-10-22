function I = loadChannels(settings, iFrame)

% Loading the sequence depends on what other kind of data need to be
% displayed.

nChannels = numel(settings.channels);
nMasks = numel(settings.maskFileNames);

% Get the channel plugins list
channelPlugins = getPlugins();

if nMasks
  BW = imread(fullfile(settings.maskPath,settings.maskFileNames{iFrame}));
end

if nChannels
  I = [];
  colors = zeros(3, 1);
  for iChannel = 1:nChannels
    channelTypeID = settings.channels{iChannel}.type;
    channelColor = settings.channels{iChannel}.color;
    
    fileName = fullfile(settings.channels{iChannel}.path,...
      settings.channels{iChannel}.fileNames{iFrame});
    
    J = channelPlugins(channelTypeID).loadFunc(fileName);
    I = cat(3, I, J);
    
    colors(iChannel) = channelColor;
  end
  
  [m n c] = size(I);
  
  if c > 1
    % RGB images are 3 x [0..1] so each channel has to be
    % normalized.
    I = single(I);
    for iChannel = 1:c
      J = I(:, :, iChannel);
      J = J - min(J(:));
      I(:, :, iChannel) = J / max(J(:));
    end
    
    if c == 2
      % Add a blank channel.
      I = cat(3, I, zeros(m, n, class(I)));
      % Set the color of the blank channel to the free color (i.e.
      % different form colors(1) and colors(2)).
      colors(3) = 6 - (colors(1) + colors(2));
    end
    
    % Make channel permutation according to channel colors
    I = I(:, :, colors);
  end
end

if nMasks
  if nChannels
    idx = cell(c,1);
    idx{1} = find(BW == 0);
    for iChannel = 2:size(I,3)
      idx{2} = idx{1} + numel(BW);
    end
    I(vertcat(idx{:})) = 0;
  else
    I = uint8(BW);
  end
  clear BW;
end
