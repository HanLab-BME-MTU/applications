function varargout = bandsDetectGUI(varargin)
% BANDSDETECTGUI MATLAB code for bandsDetectGUI.fig
%      BANDSDETECTGUI, by itself, creates a new BANDSDETECTGUI or raises the existing
%      singleton*.
%
%      H = BANDSDETECTGUI returns the handle to a new BANDSDETECTGUI or the handle to
%      the existing singleton*.
%
%      BANDSDETECTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BANDSDETECTGUI.M with the given input arguments.
%
%      BANDSDETECTGUI('Property','Value',...) creates a new BANDSDETECTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bandsDetectGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bandsDetectGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to markerSet (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bandsDetectGUI

% Last Modified by GUIDE v2.5 21-Oct-2016 12:12:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bandsDetectGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @bandsDetectGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before bandsDetectGUI is made visible.
function bandsDetectGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bandsDetectGUI (see VARARGIN)

% Choose default command line output for bandsDetectGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bandsDetectGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bandsDetectGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadImage.
function loadImage_Callback(hObject, eventdata, handles)
% hObject    handle to loadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.note, 'String', '');

handles.imagePath = loadImageFile();
image = imread(handles.imagePath);

% Convert input to single frame 2D gray scale image
if size(image,3) == 4
    image(:,:,4) = [];
else if size(image,3) == 3
        image = rgb2gray(image);
    end
end

% handles.imInput = image;
axes(handles.axes1);
imshow(image,[]);
title('Input Image');

set(handles.adjImage, 'Enable', 'on');
set(handles.markerSet, 'Enable', 'on');
guidata(hObject,handles);

% --- Executes on button press in adjImage.
function adjImage_Callback(hObject, eventdata, handles)
% hObject    handle to adjImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

imagePath = handles.imagePath;
image = imread(imagePath);
axes(handles.axes2);
imshow(image,[]);

set(handles.note, 'String', 'Please crop the region of interest');
imInput=imcrop(handles.axes2);

% Imcomplement the image if its background is white
% IM, fineIM, fIM always has black background
IM = mat2gray(imInput);
if mean(IM(:)) > 0.5
    IM=imcomplement(mat2gray(IM));
end

adjIM = imadjust(IM,[min(IM(:)),max(IM(:))], [0,1]);
[X,Y] = meshgrid(1:size(adjIM,2), 1:size(adjIM,1));
[fineX,fineY] = meshgrid(1:.2:size(adjIM,2), 1:.2:size(adjIM,1));
fineIM = interp2(X,Y,adjIM,fineX,fineY);

axes(handles.axes2);
imshow(imcomplement(fineIM),[]);
title('Cropped Image');

handles.imageIn = fineIM;

set(handles.run, 'Enable', 'on');
set(handles.note, 'String', '');

guidata(hObject,handles);



% --- Executes on button press in markerSet.
function markerSet_Callback(hObject, eventdata, handles)
% hObject    handle to markerSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fIM = handles.imageIn;
bandMap = handles.bandMap;
markerBoundary = handles.markerBoundary;
set(handles.note, 'String', '');

prompt = {'Enter the number of marker bands:', 'Enter the marker size (From top to bottom, separated by comma):'};
markerSet = inputdlg(prompt, 'Input');
markerNum = str2num(markerSet{1});
markerSize = str2num(markerSet{2});
markerPos = zeros(markerNum, 1);
markerCount = 0;
for markerY = 1:size(fIM,1)
    for markerX = 1:round(markerBoundary)
        if bandMap(markerY, markerX) == 1
            markerCount = markerCount + 1;
            markerPos(markerCount) = markerY;
        end
    end
end

if markerCount ~= markerNum
    error('Incorrect marker bands number')
end

% markerPos = markerPos';

switch markerNum
    case 7
        markerSize = [9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
    case 8
        markerSize = [18.8, 9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
    otherwise
        set(handles.note, 'String', 'Unusal number of marker bands detected');
        return
end

shortThresh = inputdlg('Enter the short telomere threshold:', 'Input');
handles.shortThresh = str2num(shortThresh{:});
handles.markerNum = markerNum;
handles.markerPos = markerPos;
handles.markerSize = markerSize;
set(handles.calc, 'Enable', 'on');
guidata(hObject,handles);



% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.note, 'String', '');
fIM = handles.imageIn;
% Lanes and bands detection

bandMap = zeros(size(fIM));
% Project all pixels intensity to x axis
intensityProfile = zeros(1,size(fIM,2));
for j = 1:size(fIM,2)
    for i = 1:size(fIM,1)
        intensityProfile(j) = intensityProfile(j) + fIM(i,j);
    end
end
% figure, plot(intensityProfile),title('Intensity Profile of filtered image')

% Should all range relavant calculations be normalized???
% Instead of pointing out the center, can you watershed out a range that
% indicates bands???

% Lanes and bands detection.
% Carefully adjust laneThresh and bandThresh
laneThresh = 0.01*(max(intensityProfile)-min(intensityProfile));
laneCenterLoc = find(~watershed(imhmin(intensityProfile,laneThresh)));
% disp(strcat(num2str(length(laneCenterLoc)), ' lane(s) detected'))

halfBandRange = max(1, round(0.3*(median(diff(laneCenterLoc))/2)));
for k = 1:length(laneCenterLoc)
    % Take average intersity value of laneCenter +/- halfBandRange
    bandProfile = zeros(size(fIM,1),1);
    laneLeftBoundary = max(0,laneCenterLoc(k)-halfBandRange);
    laneRightBoundary = min(size(fIM,2),laneCenterLoc(k)+halfBandRange);
    count = 0;
    for bandRange = laneLeftBoundary:laneRightBoundary
        bandProfile = bandProfile + fIM(:,bandRange);
        count = count + 1;
    end
    bandProfile = bandProfile/count;
    %figure,plot(laneCenter)
    %findpeaks(laneCenter)
    
    % Needs smarter threshold to eliminate the noise without hurting peaks
    % Determine the sensitivity of bands detection
    bandThresh = 0.03*(max(bandProfile)-min(bandProfile));
    bandLoc = find(~watershed(imhmin(bandProfile,bandThresh)));
    
    for m = 1:length(bandLoc)
        bandMap(bandLoc(m),laneCenterLoc(k)) = 1;
    end
    
    clear bandLoc

end

% Show preliminary detected result for manual modification
bandMapPlot(fIM, bandMap, handles);

impixelinfo;
% Marker lane is required to be on the left
set(handles.note, 'String', 'Click to define the boundary of marker lane.');
[markerBoundary, limitXX] = ginput(1);
hold on
yAxis = ylim;
plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')

% Reset marker lane in case markers are not detected??
handles.bandMap = bandMap;
handles.markerBoundary = markerBoundary;
handles.halfBandRange = halfBandRange;
set(handles.addBands, 'Enable', 'on');
set(handles.delBands, 'Enable', 'on');
set(handles.note, 'String', '');
guidata(hObject, handles);


% --- Executes on button press in addBands.
function addBands_Callback(hObject, eventdata, handles)
% hObject    handle to addBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.note, 'String', 'Click on image to add a band');
fIM = handles.imageIn;
bandMap = handles.bandMap;
markerBoundary = handles.markerBoundary;

% Manual adjustment of detected bands
[x,y] = ginput(1);
x = round(x,0);
y = round(y,0);
bandMap(y,x)=1;
bandMapPlot(fIM, bandMap, handles);
hold on
yAxis = ylim;
plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')

handles.bandMap = bandMap;
set(handles.note, 'String', '');
guidata(hObject, handles);


% --- Executes on button press in delBands.
function delBands_Callback(hObject, eventdata, handles)
% hObject    handle to delBands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.note, 'String', 'Please choose the region with bands you want to delete and then double click');
fIM = handles.imageIn;
bandMap = handles.bandMap;
markerBoundary = handles.markerBoundary;

% Manual adjustment of detected bands
[delRegion, boundaryInfo] = imcrop(handles.axes2);
for xCord = max(1, floor(boundaryInfo(1))) : min(size(fIM, 2), ceil(boundaryInfo(1)+boundaryInfo(3)))
    for yCord = max(1, floor(boundaryInfo(2))) : min(size(fIM, 1), ceil(boundaryInfo(2)+boundaryInfo(4)))
        if bandMap(yCord, xCord) == 1
            bandMap(yCord, xCord) = 0;
        end
    end
end

bandMapPlot(fIM, bandMap, handles);
hold on
yAxis = ylim;
plot(markerBoundary * ones(1,2), [yAxis(1), yAxis(2)], 'b-')

handles.bandMap = bandMap;
set(handles.note, 'String', '');
guidata(hObject, handles);


% --- Executes on button press in calc.
function calc_Callback(hObject, eventdata, handles)
% hObject    handle to calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fIM = handles.imageIn;
bandMap = handles.bandMap;
markerBoundary = handles.markerBoundary;
halfBandRange = handles.halfBandRange;
imagePath = handles.imagePath;
set(handles.note, 'String', '');

%% Band size annotation
% Enter threshold to calculate shortest telomere ratio
% shortThresh = input('Enter the threshold marker size > ');
shortThresh = handles.shortThresh;
markerNum = handles.markerNum;
markerPos = handles.markerPos;
markerSize = handles.markerSize;

% markerNum = 0;
% for markerY = 1:size(fIM,1)
%     for markerX = 1:round(markerBoundary)
%         if bandMap(markerY, markerX) == 1
%             markerNum = markerNum + 1;
%             markerPos(markerNum) = markerY;
%         end
%     end
% end
% markerPos = markerPos';
% 
% switch markerNum
%     case 7
%         markerSize = [9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
%     case 8
%         markerSize = [18.8, 9.4, 6.1, 5.4, 4.4, 3.3, 1.6, 0.8]';
%     otherwise
%         set(handles.note, 'String', 'Unusal number of marker bands detected');
%         return
% end

limitY = markerPos(markerSize == shortThresh);
% Linear regression y=ax+b for any two adjacent marker pair (x:markerPos, y:markerSize)
G(:,1) = markerPos;
G(:,2) = ones(markerNum,1);
para = struct('slope',[],'intercept',[]);
for num = 1:markerNum-1
    regPara = G(num:num+1,:)\markerSize(num:num+1);
    para(num).slope = regPara(1);
    para(num).intercept = regPara(2);
end
% para1 = G(1:2,:)\markerSize(1:2);
% para2 = G(2:3,:)\markerSize(2:3);
% para3 = G(3:4,:)\markerSize(3:4);
% para4 = G(4:5,:)\markerSize(4:5);
% para5 = G(5:6,:)\markerSize(5:6);
% para6 = G(6:7,:)\markerSize(6:7);


% Band size annotation
bandStat = struct('index', [], 'bandPos', [], 'bandSize', [], 'bandIntensity', [], 'count', []);
bandIndex = 0;
for q=1:size(fIM,2)
    for p=1:size(fIM,1)
        if bandMap(p,q)==1 && q > markerBoundary
            bandIndex = bandIndex + 1;
            bandStat(bandIndex).index = bandIndex;
            % Bands positions are recorded in cartesian coordinate
            bandStat(bandIndex).bandPos = [q, p];
            bandStat(bandIndex).count = 1;
            
            if p <= markerPos(1)
                bandStat(bandIndex).bandSize = para(1).slope*p+para(1).intercept;
            else if p > markerPos(markerNum)
                    bandStat(bandIndex).bandSize = para(markerNum-1).slope*p+para(markerNum-1).intercept;
                else
                    for num = 1:markerNum-1
                        if p > markerPos(num) && p <= markerPos(num+1)
                            bandStat(bandIndex).bandSize = para(num).slope*p+para(num).intercept;
                        end
                    end
                end
            end
            
            % Take average intensity around a band seed
            halfBandHeight = round(0.2*halfBandRange);
            bandIntensity = 0;
            bandLeftBoundary = max(0, q - halfBandRange);
            bandRightBoundary = min(size(fIM,2), q + halfBandRange);
            bandTopBoundary = max(0, p - halfBandHeight);
            bandBotBoundary = min (size(fIM,1), p + halfBandHeight);
            count = 0;
            for horiRange = bandLeftBoundary:bandRightBoundary
                for vertiRange = bandTopBoundary:bandBotBoundary
                    bandIntensity = bandIntensity + fIM(vertiRange, horiRange);
                    count = count + 1;
                end
            end
            
            bandStat(bandIndex).bandIntensity = bandIntensity/count;
        end
    end
end

for negativeBand = 1:size(bandStat,2)
    if bandStat(negativeBand).bandSize<0
        bandStat(negativeBand).bandSize=0;
    end
end

clear bandMap

%% Overlapping bands detection and counting
% Interval setup by percentile of bandsize
minBandSize = min([bandStat(:).bandSize]);
maxBandSize = max([bandStat(:).bandSize]) + 0.0001;
step = 5;
interDistance = (maxBandSize - minBandSize) / step;
bottom = minBandSize;
for i = 1:step
    top = bottom + interDistance;
    subBandStat = bandStat([bandStat(:).bandSize] >= bottom & [bandStat(:).bandSize] < top);
    % Choice between median intensity and average intensity
    medIntensity = median([subBandStat(:).bandIntensity]);
    
    % Set threshold to select highly intensed bands within the top/bottom interval
    intensityThreshold = medIntensity*1.5;
    overlappingBand = subBandStat([subBandStat(:).bandIntensity] > intensityThreshold);
    for addCount = 1:numel(overlappingBand)
        overlappingBand(addCount).count = round(overlappingBand(addCount).bandIntensity/medIntensity);
        bandStat([bandStat(:).index] == overlappingBand(addCount).index).count = overlappingBand(addCount).count;
        %         bandIndexOld = bandIndex;
        %         bandIndex = bandIndex + overlappingBand(addCount).count - 1;
        %         for indexAdd = 1:(bandIndex - bandIndexOld)
        %             bandStat(indexAdd) = bandStat([bandStat(:).index] == overlappingBand(addCount).index);
        %             bandStat(indexAdd).count = 'Overlapping Band';
        %         end
    end
    
    bottom = top;
end

%% Calculate bands statistics

% Plot identification results

axes(handles.axes2);
% close all
% figure;
imshow(imcomplement(fIM)), hold on
% set(gcf,'position',get(0,'screensize'))

countShort = 0;
countTotal = 0;

for i = 1:numel(bandStat)
    if bandStat(i).count == 1
        plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'r.')
    else if bandStat(i).count == 2
            plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'g.')
        else if bandStat(i).count > 2
                plot(bandStat(i).bandPos(1), bandStat(i).bandPos(2), 'm.')
            end
        end
    end
    
    countTotal = countTotal + bandStat(i).count;
    if bandStat(i).bandPos(2) > limitY
        countShort = countShort + bandStat(i).count;
    end
end

[pathName, fileName, ext] = fileparts(imagePath);

% Plot a threshold line
plot(1:q,limitY*ones(1,q),'b')
s1 = sprintf('Input file name: %s%s\n', fileName, ext);

avgBandSize = [bandStat(:).bandSize]*[bandStat(:).count]'/sum([bandStat(:).count]);
s2 = sprintf('The average telomere length is %.2f kb.\n', avgBandSize);

% Calculate shortest telomeres (below 1.6kb) ratio
ratio = countShort/countTotal*100;
s3 = sprintf('The ratio of shortest telomere below %.1fkb is %.2f%%.\n',shortThresh, ratio);

% Calculate the shortest 20% band size threshold
[sortBandSize, order] = sort([bandStat(:).bandSize]);
sortCount = [bandStat(order).count];
short20Point = sum(sortCount) * 0.2;
cumulativeSumCount = cumsum(sortCount);
boundaryPoint = find(cumulativeSumCount <= short20Point, 1, 'last');
if cumulativeSumCount(boundaryPoint) == short20Point
    short20Size = mean(sortBandSize([boundaryPoint, boundaryPoint+1]));
else
    short20Size = sortBandSize(boundaryPoint+1);
end

s4 = sprintf('The shortest 20%% telomere threshold is %.2f kb.\n', short20Size);

%% fIM storage???????
save([fullfile(pathName, fileName), '.mat'],'ratio', 'avgBandSize', 'bandStat', 'fIM', 'short20Size', 'imagePath')

% plot(1:q,short20Pos*ones(1,q),'m')
% hold on
% impixelinfo;
title(fileName)
g(1) = gcf;

% Size distribution histogram
axes(handles.axes1);
% figure,
h = histogram([bandStat(:).bandSize], 'BinWidth', 1);
% Normalization
bar(h.BinEdges(2:size(h.BinEdges,2))-0.5, h.Values./sum(h.Values)*100, 1, 'b')
xlabel('Telomere size (Kb)')
ylabel('Percentage of detected bands (%)')
title([fileName, ' Telomere Size Distribution'])
g(2) = gcf;

set(handles.note, 'String', {s1, s2, s3, s4}, 'FontSize', 10);
% savefig(g,fullfile(pathName, fileName))
guidata(hObject, handles);