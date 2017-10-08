function [winTrack, x, y, jumpFr] = checkWindowJump(MD, varargin)
% checkWindowJump Plot the trajectory of the 1st window in the 1st layer in
% order to find frames where window jumping happened. It plots the
% trajectory, the cell boundary at the last frame and estimated jump frames. 
%
% Usage:
%       [winTrack, x, y, jumpFr] = checkWindowJump(MD)
%
% Input:
%       MD          - a movieData object
% Output:
%       winTrack    - a matlab figure object showing the trajectory
%       x, y        - pixel coordinates of windows{1}{1} over time frame
%       jumpFr      - time frames detected when window jumping happened
%
% Option:
%       figFlag     - if 'on', matlab figures are ploted. Default is 'off'.
%
% Jungsik Noh, 2016/11/14


ip = inputParser;
ip.addParameter('figFlag', 'off');
parse(ip, varargin{:})
p = ip.Results;
figFlag = p.figFlag;

%% Input

%Make sure the movie has been windowed, and find the desired process.
iWinProc = MD.getProcessIndex('WindowingProcess');
winProc= MD.processes_{iWinProc};

iSegProc = winProc.funParams_.SegProcessIndex;   
segProc = MD.processes_{iSegProc};

MaskChannelIndex = segProc.funParams_.ChannelIndex(1);

%% last frame boundary, 1st, last window{1}{1}

currMask = segProc.loadChannelOutput(MaskChannelIndex, MD.nFrames_);

Bd = bwmorph(currMask, 'remove');

winTrack = figure('Visible', figFlag);

imagesc(Bd);
colormap(jet)
axis off

hold on

window1 = winProc.loadChannelOutput(1, MaskChannelIndex);  % (iChan, iFrame)
windowEnd = winProc.loadChannelOutput(MD.nFrames_, MaskChannelIndex);  % (iChan, iFrame)

%plotString = {'red','FaceAlpha',0};
if ~isempty(window1{1}) 
if ~isempty(window1{1}{1}) 
    plotWindows(window1{1}{1}); %, plotString)   % windows{w}{l}
end
end

if ~isempty(windowEnd{1}) 
if ~isempty(windowEnd{1}{1}) 
    plotWindows(windowEnd{1}{1}); %, plotString)   % windows{w}{l}
end
end


%% win trajectory

x = nan(MD.nFrames_, 1);
y = nan(MD.nFrames_, 1);
for fr = 1:MD.nFrames_
    winTmp = winProc.loadChannelOutput(fr, MaskChannelIndex);
    if ~isempty(winTmp{1})
        x(fr) = winTmp{1}{1}{end}(1, end); 
        y(fr) = winTmp{1}{1}{end}(2, end); 
    end
end

plot(x(~isnan(x)), y(~isnan(y)), 'y')

%% jump detection

dx = diff(x); dy = diff(y);
diffDist = sqrt(dx.^2 + dy.^2);
diffDist = diffDist(~isnan(diffDist));
%diffDistZ = nanZscore(diffDist); % Not robust

robustZ = (diffDist - median(diffDist, 'omitnan')) / mad(diffDist, 1) / 1.4826;  % see mad.m
jumpFr0 = find( abs(robustZ) > 10);
jumpFr = jumpFr0 + 1;

if numel(jumpFr) > 0
    for k = 1:numel(jumpFr)
        text(x(jumpFr(k)), y(jumpFr(k)), ['fr', num2str(jumpFr(k))], 'Color', 'red');
    end
end

title(['The 1st window trajectory. JumpFr: ', num2str(jumpFr(:)')]);


end


