function varargout = MovThresh(varargin)
% MOVTHRESH MATLAB code for MovThresh.fig
%      MOVTHRESH, by itself, creates a new MOVTHRESH or raises the existing
%      singleton*.
%
%      H = MOVTHRESH returns the handle to a new MOVTHRESH or the handle to
%      the existing singleton*.
%
%      MOVTHRESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOVTHRESH.M with the given input arguments.
%
%      MOVTHRESH('Property','Value',...) creates a new MOVTHRESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MovThresh_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MovThresh_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MovThresh

% Last Modified by GUIDE v2.5 14-Mar-2013 21:35:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MovThresh_OpeningFcn, ...
                   'gui_OutputFcn',  @MovThresh_OutputFcn, ...
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
function MovThresh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MovThresh (see VARARGIN)

% Choose default command line output for MovThresh
handles.output = hObject;


% Global variables
handles.IMAGECELL = 0; % Cell of thresholded images
handles.NUM_ELEMENTS = 0; % Number of frames in the multi-page tiff file
handles.IMAGE_NUMBER = 1; % Position of the current image in the image cell
handles.MEAN = 0; % Mean threshold level of all frames
handles.CUST_VALUES = 0;
handles.CUST_X=0;
handles.CUST_Y=0;
handles.CUST_LINE = 0;
handles.GAUSSIAN = 0;
handles.DEFAULT_GAUSSIAN = 0;
handles.LastSel = 1;
handles.REDLINE = 0; 
handles.X = 0; % Array of time points.
handles.Y = 0; % Array of actual threshold levels by frame.
handles.BLCOL='y';
handles.Algorithm = 1;
handles.ON_CURRENT_FRAME = 0;
handles.MIN_THRESH_RANGE = 0;
handles.MAX_THRESH_RANGE = 1;

handles.SMOO = 15;

% Handle listeners
handles.lismove1=addlistener(handles.movieSlider, 'Value', 'PostSet', @(src, event)switchImage(hObject, src, event));
handles.lismove2=addlistener(handles.movieSlider, 'Value', 'PostSet', @(src, event)replot(hObject, src, event));
handles.lisdraw1=addlistener(handles.drawerSlider, 'Value', 'PostSet', @(src, event)redrawMean(hObject, src, event));
handles.lisdraw5=addlistener(handles.drawerSlider, 'Value', 'PostSet', @(src, event)translateSmooth(hObject, src, event));
handles.lisdraw2=addlistener(handles.drawerSlider, 'Value', 'PostSet', @(src, event)readjustFrame(hObject, src, event));

% Update handles structure
colormap(handles.imageLoc,'gray');
set(handles.imageLoc,'Box','on','XTick',[],'YTick',[]);
set(handles.histogramLoc,'Box','on','XTick',[],'YTick',[]);

set(handles.threshInt, 'String', 'Threshold = ','ForegroundColor',[.5 .5 .5]);
set(handles.Frame,'String','Frame # ','ForegroundColor',[.5 .5 .5]); 
set(handles.curves,'ForegroundColor',[.5 .5 .5]);
set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
   
set(handles.movieSlider, 'Enable', 'off');
set(handles.drawerSlider, 'Enable', 'off');
set(handles.re_threshold, 'Enable', 'off');
set(handles.remove_single, 'Enable', 'off');
set(handles.reset, 'Enable', 'off');
set(handles.gaussian, 'Enable', 'off');
set(handles.custom, 'Enable', 'off');
set(handles.suggested, 'Enable', 'off');
set(handles.opt, 'Enable', 'off');
set(handles.save_as, 'Enable', 'off');
set(handles.threshold_range, 'Enable', 'off');

set(handles.boost,'Enable','off');
set(handles.edit_smoo,'Enable','off');


guidata(hObject, handles);

function varargout = MovThresh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% -------------------------------------------------------------------------



% --- Function for Listeners ----------------------------------------------
function switchImage(hObject, src, event)
% --- Switches between images in the cell of thresholded images by 
% listening to the movie slider. The sequence of images is dependent upon
% the last thresholding of the entire movie.
    
%disp('switch');

    handles = guidata(hObject);
    
    sliderNumber = round(event.NewValue);
    
    handles.IMAGE_NUMBER = sliderNumber;
    
    axes(handles.imageLoc);
    xlim=get(handles.imageLoc,'XLim');
    ylim=get(handles.imageLoc,'YLim');
    cla;
    hold on;
    imagesc(handles.IMAGECELL{handles.IMAGE_NUMBER}, 'Parent', handles.imageLoc); 
    axis image;
    axis ij;
    set(handles.imageLoc,'XLim',xlim);
    set(handles.imageLoc,'YLim',ylim);
    axis off;
    bound=handles.BOUNDARIES{handles.IMAGE_NUMBER};
    by=bound(:,1); bx=bound(:,2);
    plot(bx,by,handles.BLCOL,'LineWidth',2);
    %plot(bx,by,'m','LineWidth',2);
    
    if(get(handles.suggested, 'Value') == 1)
        thresh = handles.Y(handles.IMAGE_NUMBER);
        
    elseif(get(handles.gaussian, 'Value') == 1)
        thresh = handles.GAUSSIAN(handles.IMAGE_NUMBER);
        
    elseif(get(handles.custom, 'Value') == 1)
        thresh = handles.CUST_VALUES(handles.IMAGE_NUMBER);
    end
    
    str=sprintf('Threshold = %.2f',thresh);
    set(handles.threshInt, 'String', str);
    set(handles.Frame,'String',['Frame # ' num2str(handles.IMAGE_NUMBER)]);
    
    guidata(hObject, handles);
  
function replot(hObject, src, event)
% --- Listens to the movie slider and replots the red line on the line graph of 
% optimal intensities for every time unit.

    handles = guidata(hObject);
                 
    handles.ON_CURRENT_FRAME = 0;
    
    value = round(event.NewValue);
    
    if ishandle(handles.REDLINE) 
        if handles.REDLINE~=0
            delete(handles.REDLINE);
        end
    end
    handles.REDLINE = line([value, value], get(handles.histogramLoc, 'YLim'), 'Color', 'r', 'Parent', ...
        handles.histogramLoc, 'EraseMode', 'normal'); 
    
    if(strcmp(get(handles.drawerSlider, 'Enable'), 'on'))
        if (get(handles.custom, 'Value') == 1)
            set(handles.drawerSlider, 'Value', handles.CUST_VALUES(value));
        elseif (get(handles.gaussian, 'Value') == 1)
            set(handles.drawerSlider, 'Value', handles.gaussianSavedSliderValue);
        end
    end
    
    guidata(hObject, handles);
 
function redrawMean(hObject, src, event) 
% -- Redraws the custom curve on the intensity plot when the user changes 
% the threshold point of one frame.   


    handles = guidata(hObject);
    switchPoint = handles.IMAGE_NUMBER;
    newIntensity = event.NewValue;
    
    if (get(handles.custom, 'Value') ~= 1)
        return;
    end
    
    if(newIntensity ~= handles.CUST_VALUES(switchPoint)) 
        handles.CUST_VALUES(switchPoint) = newIntensity;
        
        if switchPoint > 1 && switchPoint < handles.NUM_ELEMENTS
            i=find(handles.CUST_X>switchPoint,1,'first');
            j=find(handles.CUST_X<switchPoint,1,'last');
            i=handles.CUST_X(i);
            j=handles.CUST_X(j);
        elseif switchPoint == 1
            i=find(handles.CUST_X>switchPoint,1,'first');
            i=handles.CUST_X(i);
            j=1;
        else
            i=handles.NUM_ELEMENTS;
            j=find(handles.CUST_X<switchPoint,1,'last');
            j=handles.CUST_X(j);
        end        
        
        slope1 = (newIntensity - handles.CUST_VALUES(j))/(switchPoint - j);
        intercept1 = newIntensity - (slope1.*switchPoint);

        slope2 = (handles.CUST_VALUES(i) - newIntensity)/(i - switchPoint);
        intercept2 = newIntensity - (slope2.*switchPoint);

        if(switchPoint ~= 1)
            for k = j+1:switchPoint-1
                handles.CUST_VALUES(k) = slope1.*k + intercept1;
            end
        end
        
        if(switchPoint ~= handles.NUM_ELEMENTS)
            for k = switchPoint+1:i-1
               handles.CUST_VALUES(k) = slope2.*k + intercept2;
            end
        end
        
        q=find(handles.CUST_X==switchPoint,1);
        if isempty(q)
            tmp=[handles.CUST_Y newIntensity];
            [handles.CUST_X,ind]=sort([handles.CUST_X switchPoint]);
            handles.CUST_Y=tmp(ind);
        else
            handles.CUST_Y(q)=newIntensity;
        end
                        
        axes(handles.histogramLoc);
        if ~ishandle(handles.CUST_LINE) || handles.CUST_LINE == 0
            return;
        end
        delete(handles.CUST_LINE);
        handles.CUST_LINE = line('XData', handles.CUST_X, 'YData', handles.CUST_Y, 'Color', 'k','Marker','o','MarkerFaceColor','r','MarkerSize',4); 
        str=sprintf('Threshold = %.2f',handles.CUST_VALUES(handles.IMAGE_NUMBER));
        set(handles.threshInt, 'String', str);
    end

    guidata(hObject, handles);

function translateSmooth(hObject, src, event)

    handles = guidata(hObject);
    newValue = event.NewValue;
    
    if (get(handles.gaussian, 'Value') ~= 1)
        return;
    end
    
    if (newValue ~= handles.gaussianSavedSliderValue)
        handles.gaussianSavedSliderValue = newValue;
        X = handles.X;
        Y = handles.DEFAULT_GAUSSIAN-(handles.MEAN-newValue);
        handles.GAUSSIAN = Y;

        if ishandle(handles.CUST_LINE) && handles.CUST_LINE ~= 0
            delete(handles.CUST_LINE);
        end
        
        axes(handles.histogramLoc);
        handles.CUST_LINE = line('Xdata', X, 'YData', Y, 'Color', 'r', 'EraseMode', 'normal','Marker','.');
        str=sprintf('Threshold = %.2f',Y(handles.IMAGE_NUMBER));
        set(handles.threshInt, 'String', str);
    end
    
    guidata(hObject, handles);
    
    
function readjustFrame(hObject, src, event)
% --- If the drawer slider is moved directly after thresholding a single 
% frame, the threshold boundary on the current frame will change according 
% to the threshold value on the intensity plot.
    
    handles = guidata(hObject);
    
    if handles.ON_CURRENT_FRAME == 1  
        value = event.NewValue;

        upper = get(handles.drawerSlider, 'Max');
        lower = get(handles.drawerSlider, 'Min');

        axes(handles.imageLoc);
        xlim=get(handles.imageLoc,'XLim');
        ylim=get(handles.imageLoc,'YLim');
        %cla(handles.imageLoc, 'reset'); 
        cla;
        hold on;

        image = handles.IMAGECELL{handles.IMAGE_NUMBER};
        %image = imread(handles.FILENAME, handles.IMAGE_NUMBER);
        imagesc(image);
        index = ((value-lower)/(upper-lower)).*100 + 1;
        bound=handles.CURRENT_BOUND{round(index)};
        by=bound(:,1); bx=bound(:,2);
        plot(bx,by,handles.BLCOL,'LineWidth',2);                
        handles.BOUNDARIES{handles.IMAGE_NUMBER} = bound;
        
        axis image;
        axis ij;
        set(handles.imageLoc,'XLim',xlim);
        set(handles.imageLoc,'YLim',ylim);
        axis off;        
        
        thresh= lower + (upper-lower) * (round(index)-1)/100;
        str=sprintf('Threshold = %.2f',thresh);
        set(handles.threshInt, 'String', str);
    
    end
    
    guidata(hObject, handles);
% -------------------------------------------------------------------------



% --- Additional functions ------------------------------------------------
function import(hObject,boost)
% --- Imports a multi-page tiff image from a given directory and thresholds it by frame. 
    jFigPeer = get(handle(gcf),'JavaFrame');
    try
        jWindow = jFigPeer.fFigureClient.getWindow;
    catch
        jWindow = jFigPeer.fHG1Client.getWindow;
    end
    set(handle(jWindow),'Enabled',false);
    try

        handles = guidata(hObject);
        
        %try                
        %    delete(handles.lisdraw1);          
        %end
        %try                
        %    delete(handles.lisdraw2);          
        %end
        
        axes(handles.histogramLoc);
        cla;

        axes(handles.imageLoc);
        cla reset; %axis fill;
        %set(handles.imageLoc,'Position',[20 62 500 430]);
        set(handles.imageLoc,'Box','on','XTick',[],'YTick',[],'Visible','on');

        set(handles.threshInt, 'String', 'Threshold = ','ForegroundColor',[.5 .5 .5]);
        set(handles.Frame,'String','Frame # ','ForegroundColor',[.5 .5 .5]); 
        set(handles.curves,'ForegroundColor',[.5 .5 .5]);

        set(handles.movieSlider, 'Enable', 'off');
        set(handles.drawerSlider, 'Enable', 'off');
        set(handles.re_threshold, 'Enable', 'off');
        set(handles.remove_single, 'Enable', 'off');
        set(handles.reset, 'Enable', 'off');
        set(handles.gaussian, 'Enable', 'off');
        set(handles.custom, 'Enable', 'off');
        set(handles.suggested, 'Enable', 'off');
        set(handles.opt, 'Enable', 'off');
        set(handles.save_as, 'Enable', 'off');
        set(handles.threshold_range, 'Enable', 'off');
        set(handles.boost,'Enable','off');
        set(handles.suggested, 'Value', 1);
        
        set(handles.lismove1, 'Enable', 'off');
        set(handles.lismove2, 'Enable', 'off');
        set(handles.lisdraw1, 'Enable', 'off');
        set(handles.lisdraw2, 'Enable', 'off');
        

        [filename,pathname,filterindex] = uigetfile('*.tif'); 

        if filename~=0
            
            handles.ON_CURRENT_FRAME = 0;
            handles.IMAGE_NUMBER = 1;

            info = imfinfo([pathname filename]);
            handles.NUM_ELEMENTS = numel(info);   
                        
            if handles.NUM_ELEMENTS>1

                set(handles.movieSlider, 'SliderStep',[1/(handles.NUM_ELEMENTS-1) 5/(handles.NUM_ELEMENTS-1)]);
                set(handles.movieSlider, 'Min', 1);
                set(handles.movieSlider, 'Max', handles.NUM_ELEMENTS);
                set(handles.histogramLoc, 'XLim', [1 handles.NUM_ELEMENTS]);   

                if handles.SMOO >= handles.NUM_ELEMENTS
                    handles.SMOO = handles.NUM_ELEMENTS-1;
                end 

                ImageCell = cell(1,handles.NUM_ELEMENTS);

                x = 1:handles.NUM_ELEMENTS;
                y = zeros(1,handles.NUM_ELEMENTS);
                b = cell(1,handles.NUM_ELEMENTS);
                
                timage = imread([pathname filename], 1);
                if length(size(timage))==3
                    timage=timage(:,:,1);
                end
                if boost
                    timage=EnhanceGrid(timage);                   
                end
                set(handles.imageLoc,'XLim',[0.5 size(timage,2)+0.5]);
                set(handles.imageLoc,'YLim',[0.5 size(timage,1)+0.5]);
                

                for k = 1:handles.NUM_ELEMENTS
                   axes(handles.progress1);
                   cla;                 
                   bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
                   by=[0 0 1 1 0];
                   fill(bx,by,[.7 .7 .7]);
                   set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
                   pause(0.001); 
                   %set(handles.Frame,'String',['Frame # ' num2str(k)]); 
                   image = imread([pathname filename], k);
                   if length(size(image))==3
                        image=image(:,:,1);
                   end
                   image = cast(image, 'double');
                   if boost
                       image=EnhanceGrid(image);                   
                   end
                   ImageCell{k}=image;
                   [y(k), b{k}] = threshold(hObject,image);             
                end
                cla(handles.progress1,'reset');   
                set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');

                handles.IMAGECELL = ImageCell;
                handles.X = x;
                handles.Y = y;
                handles.BOUNDARIES = b;
                handles.IMAGE_NUMBER = 1;

                gaussianY = GFilterA(y, handles.SMOO); 
                handles.GAUSSIAN = gaussianY;
                handles.DEFAULT_GAUSSIAN = gaussianY;

                handles.MEAN = mean(y);
                handles.CUST_VALUES = repmat(handles.MEAN, 1, handles.NUM_ELEMENTS);
                handles.CUST_X=[1 handles.NUM_ELEMENTS];
                handles.CUST_Y=[handles.MEAN handles.MEAN];

                axes(handles.histogramLoc);
                hold on;
                plot(x, y, '.-'); 
                range=max(y)-min(y);
                if range>0
                    set(handles.histogramLoc, 'YLim', [floor(min(y)-range/3) ceil(max(y)+range/3)]);
                else                
                    set(handles.histogramLoc, 'YLim', [floor(max(y)-1) ceil(max(y)+1)]);
                end       

                if ~ishandle(handles.REDLINE) || handles.REDLINE==0
                    handles.REDLINE = line([1, 1], get(handles.histogramLoc, 'YLim'), 'Color', 'r', 'Parent', ...
                        handles.histogramLoc, 'EraseMode', 'normal'); 
                end

                axes(handles.imageLoc);
                cla;
                hold on;
                imagesc(ImageCell{1});
                bound=b{1};
                by=bound(:,1); bx=bound(:,2);
                plot(bx,by,handles.BLCOL,'LineWidth',2);
                axis image; 
                axis ij;
                axis off;
                zoom reset;
                
                set(handles.movieSlider, 'Value', 1);
                                
                if range>0
                    minv = floor(min(y)-range/3);
                    maxv = ceil(max(y)+range/3);
                else                
                    minv = floor(max(y)-1);
                    maxv = ceil(max(y)+1);
                end 
                set(handles.drawerSlider, 'Min', minv, 'Max', maxv);
                set(handles.min_thresh_range, 'String', num2str(minv));
                set(handles.max_thresh_range, 'String', num2str(maxv));
                handles.MIN_THRESH_RANGE = minv;
                handles.MAX_THRESH_RANGE = maxv;
                
                str=sprintf('Threshold = %.2f',y(1));
                set(handles.threshInt, 'String', str,'ForegroundColor','k');
                set(handles.Frame,'String','Frame # 1','ForegroundColor','k');

                %guidata(hObject, handles);
                set(handles.drawerSlider, 'Value', handles.MEAN);
                
                % save the slider value for gaussian smooth curve
                handles.gaussianSavedSliderValue = handles.MEAN;
                
                %handles.lisdraw1=addlistener(handles.drawerSlider, 'Value', 'PostSet', @(src, event)redrawMean(hObject, src, event));
                %handles.lisdraw2=addlistener(handles.drawerSlider, 'Value', 'PostSet', @(src, event)readjustFrame(hObject, src, event));

                set(handles.movieSlider, 'Enable', 'on');
                set(handles.re_threshold, 'Enable', 'on');
                set(handles.gaussian, 'Enable', 'on');
                set(handles.custom, 'Enable', 'on');
                set(handles.suggested, 'Enable', 'on');
                set(handles.opt, 'Enable', 'on');
                set(handles.save_as, 'Enable', 'on');
                set(handles.edit_smoo, 'Enable', 'off', 'String', ' ');
                
                set(handles.lismove1, 'Enable', 'on');
                set(handles.lismove2, 'Enable', 'on');
                set(handles.lisdraw1, 'Enable', 'on');
                set(handles.lisdraw2, 'Enable', 'on');
                
            elseif handles.NUM_ELEMENTS==1
                
                set(handles.Frame,'String',['Frame # ' num2str(1)],'ForegroundColor','k'); 
                image = imread([pathname filename], 1);
                if length(size(image))==3
                    image=image(:,:,1);
                end
                image = cast(image, 'double');
                b=cell(1,1);
                if boost
                    image=EnhanceGrid(image);                   
                end
                handles.IMAGECELL=cell(1,1);
                handles.IMAGECELL{1}=image;
                
                [y, b{1}] = threshold(hObject,image);                
                                
                str=sprintf('Threshold = %.2f',y);
                set(handles.threshInt, 'String', str,'ForegroundColor','k');
                set(handles.Frame,'String','Frame # 1');
                
                handles.X = 1;
                handles.Y = y;
                handles.BOUNDARIES = b;
                
                axes(handles.imageLoc);
                cla;                
                hold on;
                imagesc(image);
                bound=b{1};
                by=bound(:,1); bx=bound(:,2);                 
                plot(bx,by,handles.BLCOL,'LineWidth',2);                                
                axis image; 
                axis ij;
                axis off; 
                zoom reset;
                                
                minv=max(floor(0.9*max(y)),0);
                maxv=max(ceil(1.1*max(y)),1);
                
                set(handles.drawerSlider, 'Min', minv, 'Max', maxv);
                set(handles.drawerSlider, 'Value', y);    
                set(handles.min_thresh_range, 'String', num2str(minv));
                set(handles.max_thresh_range, 'String', num2str(maxv));
                handles.MIN_THRESH_RANGE = minv;
                handles.MAX_THRESH_RANGE = maxv;
                                
                set(handles.opt, 'Enable', 'on');
                set(handles.save_as, 'Enable', 'on');
                set(handles.threshold_range, 'Enable', 'on');
                set(handles.re_threshold, 'Enable', 'on');
                set(handles.edit_smoo, 'Enable', 'off', 'String', ' ');
                
                set(handles.lismove1, 'Enable', 'on');
                %set(handles.lismove2, 'Enable', 'on');
                %set(handles.lisdraw1, 'Enable', 'on');
                set(handles.lisdraw2, 'Enable', 'on');
                
            end 
            set(handles.boost,'Enable','on');    
        end

        
        
        handles.FILENAME = [pathname filename];
        guidata(hObject,handles);   
        
    catch
        errordlg('Import failed. Please check the file format.');
    end
    
    set(handle(jWindow),'Enabled',true);

function [Thresh, bound] = threshold(hObject, image)
% --- Thresholds a given image with the optimal threshold level, then displays it on the GUI.

    handles = guidata(hObject);
%----------------------< Finding a threshold >---------------------------

    if handles.Algorithm == 1
        SmooWind=10; Gsize=[5 5]; Gsigm=1;                 % smoothing parameters
        h=fspecial('gaussian', Gsize, Gsigm);
        T=imfilter(image(image>0), h,'replicate');
        [Thresh,Status,IntMax,Xhist,Yhist]=AutoHist(T,SmooWind);      %!!!!!!!!!!!
    elseif handles.Algorithm == 2     
        level=graythresh(uint8(image));
        Thresh=level*255;       
    end           
       
    %----------------------< Thresholding the image >------------------------

    ThreshImage=image-Thresh;
    ThreshImage(ThreshImage<0)=0;
    BinaryImage=imfill(logical(ThreshImage),'holes');
    if max(BinaryImage(:))==min(BinaryImage(:))
        BinaryImage=ones(size(image));
    end    
    
    %----------------------< Boundaries and Areas >---------------------------

    B=bwboundaries(BinaryImage);
    A=regionprops(logical(BinaryImage),'Area');
    Area=[A.Area];
    [v,ind]=max(Area);
    bound=B{ind};

    %----------------------< Plotting >---------------------------------------  
    %%{
    axes(handles.imageLoc);
    cla;
    hold on;
    imagesc(image);     
    by=bound(:,1); bx=bound(:,2);
    plot(bx,by,handles.BLCOL,'LineWidth',2);
    axis image;
    axis ij;
    axis off;
    %}
    
    guidata(hObject, handles);

function bound=thresholdByValue(hObject, image, intensity)
% --- Thresholds a given image with the specified intensity value, then displays it on the GUI.

    handles = guidata(hObject);
    %----------------------< Thresholding the image >------------------------

    ThreshImage=image-intensity;
    ThreshImage(ThreshImage<0)=0;
    BinaryImage=imfill(logical(ThreshImage),'holes');
    if max(BinaryImage(:))==min(BinaryImage(:))
        BinaryImage=ones(size(image));
    end
    %----------------------< Boundaries and Areas >---------------------------

    B=bwboundaries(BinaryImage);
    A=regionprops(logical(BinaryImage),'Area');
    Area=[A.Area];
    [v,ind]=max(Area);
    bound=B{ind};

    %----------------------< Plotting >---------------------------------------
    %%{
    axes(handles.imageLoc);
    xlim=get(handles.imageLoc,'XLim');
    ylim=get(handles.imageLoc,'YLim');
    cla;
    hold on;
    imagesc(image);
    axis image;
    axis ij;
    set(handles.imageLoc,'XLim',xlim);
    set(handles.imageLoc,'YLim',ylim);
    axis off;    
    by=bound(:,1); bx=bound(:,2);
    plot(bx,by,handles.BLCOL,'LineWidth',2);
    %}
    
    guidata(hObject, handles);
      
function thresholdBySuggested(hObject)
% --- Thresholds all frames in the multipage tiff image with the mean 
% threshold level of all the frames.

    handles = guidata(hObject);
    
    if(get(handles.custom, 'Value') == 1)
        suggested = handles.CUST_VALUES;
        handles.LastSel=3;   
        set(handles.reset, 'Enable', 'off');
        set(handles.remove_single, 'Enable', 'off');
        set(handles.drawerSlider, 'Enable', 'off');
        set(handles.threshold_range, 'Enable', 'off');
    elseif(get(handles.gaussian, 'Value') == 1) 
        suggested = handles.GAUSSIAN;
        handles.LastSel=2;
    elseif(get(handles.suggested, 'Value') == 1)
        suggested = handles.Y;
        handles.LastSel=1;
    end
    
    b = cell(1,handles.NUM_ELEMENTS);
    
    for k = 1:handles.NUM_ELEMENTS
       axes(handles.progress1);
       cla;     
       bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
       by=[0 0 1 1 0];
       fill(bx,by,[.7 .7 .7]);
       set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
       pause(0.001); 
       %set(handles.Frame,'String',['Frame # ' num2str(k)]);
       image = handles.IMAGECELL{k};
       b{k}=thresholdByValue(hObject, image, suggested(k));             
    end
    cla(handles.progress1,'reset');    
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
    
    handles.BOUNDARIES = b;
    
    %{
    xlim=get(handles.imageLoc,'XLim');
    ylim=get(handles.imageLoc,'YLim');
    cla;
    hold on;
    imagesc(handles.IMAGECELL{handles.IMAGE_NUMBER}, 'Parent', handles.imageLoc); 
    axis image;
    axis ij;
    set(handles.imageLoc,'XLim',xlim);
    set(handles.imageLoc,'YLim',ylim);
    axis off;
    %}
    
    axes(handles.imageLoc);
    xlim=get(handles.imageLoc,'XLim');
    ylim=get(handles.imageLoc,'YLim');
    cla;
    hold on;
    imagesc(handles.IMAGECELL{handles.IMAGE_NUMBER});   %  !!!!
    axis image;
    axis ij;
    set(handles.imageLoc,'XLim',xlim);
    set(handles.imageLoc,'YLim',ylim);
    axis off;
    bound=b{handles.IMAGE_NUMBER};                     % !!!!
    by=bound(:,1); bx=bound(:,2);
    plot(bx,by,handles.BLCOL,'LineWidth',2);
        
    
    if(get(handles.custom, 'Value') == 1)   
        set(handles.reset, 'Enable', 'on');
        set(handles.remove_single, 'Enable', 'on');
        set(handles.drawerSlider, 'Enable', 'on');
        set(handles.threshold_range, 'Enable', 'on');
    end
    
    guidata(hObject, handles);
% -------------------------------------------------------------------------



% --- Sliders -------------------------------------------------------------
function movieSlider_Callback(hObject, eventdata, handles)

function movieSlider_CreateFcn(hObject, eventdata, handles)

    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
 
function drawerSlider_Callback(hObject, eventdata, handles)

function drawerSlider_CreateFcn(hObject, eventdata, handles)
    
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
% -------------------------------------------------------------------------



% --- Buttons -------------------------------------------------------------
function remove_single_Callback(hObject, eventdata, handles)
% --- Executes on button press in remove_single. If the red line on the
% intensity plot is hovering over a pivot point, the pivot will be removed
% and a straight line will be plotted between the two nearest pivots in
% either direction.
    if(handles.IMAGE_NUMBER ~= 1 && handles.IMAGE_NUMBER ~= handles.NUM_ELEMENTS)
        
        q=find(handles.CUST_X==handles.IMAGE_NUMBER,1);
        if ~isempty(q)
            
            Xa=handles.CUST_X(q-1); Ya=handles.CUST_Y(q-1);
            Xb=handles.CUST_X(q+1); Yb=handles.CUST_Y(q+1);
            V=Ya+(handles.IMAGE_NUMBER-Xa)*(Yb-Ya)/(Xb-Xa);
                        
            handles.CUST_VALUES(Xa:Xb) = linspace(Ya,Yb,Xb-Xa+1);
            
            handles.CUST_X=[handles.CUST_X(1:(q-1)) handles.CUST_X((q+1):end)];
            handles.CUST_Y=[handles.CUST_Y(1:(q-1)) handles.CUST_Y((q+1):end)];
            
            
            axes(handles.histogramLoc);
            if ishandle(handles.CUST_LINE)
                if(handles.CUST_LINE ~= 0)
                    delete(handles.CUST_LINE);
                end
            end
            handles.CUST_LINE = line('Xdata', handles.CUST_X, 'YData', handles.CUST_Y, 'Color', 'k', 'EraseMode', 'normal','Marker','o','MarkerFaceColor','r','MarkerSize',4);  
            str=sprintf('Threshold = %.2f',V);
            set(handles.threshInt, 'String', str);
            
            guidata(hObject, handles);
            set(handles.drawerSlider, 'Value', handles.CUST_VALUES(handles.IMAGE_NUMBER));
            
        end
    
    end
        
function reset_Callback(hObject, eventdata, handles)
% --- Executes on button press in reset. Resets the custom curve on the 
% intensity plot to a horizontal line where the intensity at all times is 
% the mean threshold point.

    for k = 1:handles.NUM_ELEMENTS
        handles.CUST_VALUES(k) = handles.MEAN;
    end    
    %---DT---
    handles.CUST_X=[1 handles.NUM_ELEMENTS];
    handles.CUST_Y=[handles.MEAN handles.MEAN];
    %--------
    
    axes(handles.histogramLoc);
    if ishandle(handles.CUST_LINE)
        if(handles.CUST_LINE ~= 0)
            delete(handles.CUST_LINE);
        end
    end
    
    handles.CUST_LINE = line('XData', handles.CUST_X, 'YData', handles.CUST_Y, 'Color','k','Marker','o','MarkerFaceColor','r','MarkerSize',4);
    
    guidata(hObject, handles);
    set(handles.drawerSlider, 'Value', handles.CUST_VALUES(handles.IMAGE_NUMBER));    

function re_threshold_Callback(hObject, eventdata, handles)
% --- Executes on button press in re_threshold. Re-thresholds the entire
% sequence of images based on the curve on the intensity plot
% (user-specified).
    
    jFigPeer = get(handle(gcf),'JavaFrame');
    try
        jWindow = jFigPeer.fFigureClient.getWindow;
    catch
        jWindow = jFigPeer.fHG1Client.getWindow;
    end
    set(handle(jWindow),'Enabled',false);
    try
        
        if handles.NUM_ELEMENTS>1

            set(handles.movieSlider, 'Enable', 'off');
            set(handles.re_threshold, 'Enable', 'off');

            thresholdBySuggested(hObject);

            set(handles.movieSlider, 'Enable', 'on');
            set(handles.re_threshold, 'Enable', 'on');

            %set(handles.movieSlider, 'Value', 1);
        
        elseif handles.NUM_ELEMENTS==1
            %handles.LastSel=1;           
            handles.Y=get(handles.drawerSlider, 'Value');
            guidata(hObject, handles);
        end
            
            
    end
    set(handle(jWindow),'Enabled',true);

function threshold_range_Callback(hObject, eventdata, handles)
% --- Executes on button press in threshold_range. Thresholds a single
% frame within the intensity range specified by the user.
    
    jFigPeer = get(handle(gcf),'JavaFrame');
    try
        jWindow = jFigPeer.fFigureClient.getWindow;
    catch
        jWindow = jFigPeer.fHG1Client.getWindow;
    end
    set(handle(jWindow),'Enabled',false);
    try
        if handles.NUM_ELEMENTS>1
            
            set(handles.movieSlider, 'Enable', 'off');
            set(handles.drawerSlider, 'Enable', 'off');
            set(handles.re_threshold, 'Enable', 'off');
            set(handles.remove_single, 'Enable', 'off');
            set(handles.reset, 'Enable', 'off');
            set(handles.gaussian, 'Enable', 'off');
            set(handles.custom, 'Enable', 'off');
            set(handles.suggested, 'Enable', 'off');
            set(handles.threshold_range, 'Enable', 'off');

            upper = get(handles.drawerSlider, 'Max');
            lower = get(handles.drawerSlider, 'Min');
            sliderStep = get(handles.drawerSlider, 'SliderStep');

            image = handles.IMAGECELL{handles.IMAGE_NUMBER};

            handles.CURRENT_BOUND = cell(1, 101);

            for k = 1:101 
                axes(handles.progress1);
                cla;     
                bx=[0 k k 0 0]/101;
                by=[0 0 1 1 0];
                fill(bx,by,[.7 .7 .7]);
                set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
                pause(0.001);     
                
                threshValue = lower + (upper-lower) .* (k-1) .* sliderStep(1);
                handles.CURRENT_BOUND{k} = thresholdByValue(hObject, image, threshValue);               
                %str=sprintf('Threshold = %.2f',threshValue);
                %set(handles.threshInt, 'String', str);
            end
            cla(handles.progress1,'reset');    
            set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');

            axes(handles.imageLoc);
            xlim=get(handles.imageLoc,'XLim');
            ylim=get(handles.imageLoc,'YLim');
            cla;
            hold on;
            imagesc(image);
            %bound = handles.BOUNDARIES{handles.IMAGE_NUMBER};    
            bound = thresholdByValue(hObject, image, get(handles.drawerSlider,'Value'));   % !!!!!!!!!!
            by=bound(:,1); bx=bound(:,2);
            plot(bx,by,handles.BLCOL,'LineWidth',2);
            axis image;    
            axis ij;
            set(handles.imageLoc,'XLim',xlim);
            set(handles.imageLoc,'YLim',ylim);
            axis off;

            if(get(handles.suggested, 'Value') == 1)
                thresh = handles.Y(handles.IMAGE_NUMBER);

            elseif(get(handles.gaussian, 'Value') == 1)
                thresh = handles.GAUSSIAN(handles.IMAGE_NUMBER);

            elseif(get(handles.custom, 'Value') == 1)
                thresh = handles.CUST_VALUES(handles.IMAGE_NUMBER);
            end

            str=sprintf('Threshold = %.2f',thresh);
            set(handles.threshInt, 'String', str);

            set(handles.movieSlider, 'Enable', 'on');
            set(handles.drawerSlider, 'Enable', 'on');
            set(handles.re_threshold, 'Enable', 'on');
            set(handles.remove_single, 'Enable', 'on');
            set(handles.reset, 'Enable', 'on');
            set(handles.gaussian, 'Enable', 'on');
            set(handles.custom, 'Enable', 'on');
            set(handles.suggested, 'Enable', 'on');
            set(handles.threshold_range, 'Enable', 'on');

            handles.ON_CURRENT_FRAME = 1;
        elseif handles.NUM_ELEMENTS==1
            
            upper = get(handles.drawerSlider, 'Max');
            lower = get(handles.drawerSlider, 'Min');
            sliderStep = get(handles.drawerSlider, 'SliderStep');

            image = handles.IMAGECELL{1};

            handles.CURRENT_BOUND = cell(1, 101);
            handles.ON_CURRENT_FRAME = 1;

            for k = 1:101
                axes(handles.progress1);
                cla;     
                bx=[0 k k 0 0]/101;
                by=[0 0 1 1 0];
                fill(bx,by,[.7 .7 .7]);
                set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
                pause(0.001);  
                
                threshValue = lower + (upper-lower) .* (k-1) .* sliderStep(1);
                handles.CURRENT_BOUND{k} = thresholdByValue(hObject, image, threshValue);               
                %str=sprintf('Threshold = %.2f',threshValue);
                %set(handles.threshInt, 'String', str);
            end
            cla(handles.progress1,'reset');    
            set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
            
            axes(handles.imageLoc);
            xlim=get(handles.imageLoc,'XLim');
            ylim=get(handles.imageLoc,'YLim');
            cla;
            hold on;
            imagesc(image);
            bound = handles.BOUNDARIES{1};
            by=bound(:,1); bx=bound(:,2);
            plot(bx,by,handles.BLCOL,'LineWidth',2);
            axis image;    
            axis ij;
            set(handles.imageLoc,'XLim',xlim);
            set(handles.imageLoc,'YLim',ylim);
            axis off;
            
            set(handles.drawerSlider, 'Enable', 'on');
        end

        guidata(hObject, handles);
    end
    axes(handles.histogramLoc);
    set(handle(jWindow),'Enabled',true);
% -------------------------------------------------------------------------    



% --- Selections ----------------------------------------------------------
function curves_SelectionChangeFcn(hObject, eventdata, handles)
% --- Executes when selected object is changed in suggested. Allows the
% user to change the curve fit on the intensity plot. 
    axes(handles.histogramLoc);
    
    switch get(eventdata.NewValue, 'Tag')
        case 'gaussian'
            Y = GFilterA(handles.Y, handles.SMOO);
            X = handles.X;
            handles.DEFAULT_GAUSSIAN = Y; 
            
            if ishandle(handles.CUST_LINE)
                if(handles.CUST_LINE ~= 0) 
                    delete(handles.CUST_LINE);
                end
            end
            % restore
            set(handles.drawerSlider, 'Value', handles.gaussianSavedSliderValue);
            Y = Y-(handles.MEAN-handles.gaussianSavedSliderValue);
            handles.GAUSSIAN = Y;
            handles.CUST_LINE = line('Xdata', X, 'YData', Y, 'Color', 'r', 'EraseMode', 'normal','Marker','.', 'Parent', handles.histogramLoc); 
            str=sprintf('Threshold = %.2f',handles.GAUSSIAN(handles.IMAGE_NUMBER));
            set(handles.threshInt, 'String', str);
            set(handles.remove_single, 'Enable', 'off');
            set(handles.drawerSlider, 'Enable', 'on');
            set(handles.reset, 'Enable', 'off');
            set(handles.threshold_range, 'Enable', 'off');
            
            set(handles.smoo_text,'ForegroundColor','k');
            set(handles.edit_smoo,'Enable','on','String',num2str(handles.SMOO),'Value',handles.SMOO);
            
        case 'custom'
            if ishandle(handles.CUST_LINE)
                if(handles.CUST_LINE ~= 0) 
                    delete(handles.CUST_LINE);
                end
            end
            
            %restore
            %set(handles.drawerSlider, 'Value', handles.customSavedSliderValue);
            set(handles.drawerSlider, 'Value', handles.CUST_VALUES(handles.IMAGE_NUMBER));
            %replot();
            handles.CUST_LINE = line('Xdata', handles.CUST_X, 'YData', handles.CUST_Y, 'Color', 'k', 'EraseMode', 'normal','Marker','o','MarkerFaceColor','r','MarkerSize',4, 'Parent',handles.histogramLoc); 
            str=sprintf('Threshold = %.2f',handles.CUST_VALUES(handles.IMAGE_NUMBER));
            set(handles.threshInt, 'String', str);
            set(handles.drawerSlider, 'Enable', 'on');
            set(handles.remove_single, 'Enable', 'on');
            set(handles.reset, 'Enable', 'on');
            set(handles.threshold_range, 'Enable', 'on');
            
            set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
            set(handles.edit_smoo,'Enable','off','String','  ');
            
        case 'suggested'
            if ishandle(handles.CUST_LINE)
                if(handles.CUST_LINE ~= 0) 
                    delete(handles.CUST_LINE);
                    handles.CUST_LINE = 0;
                end
            end
            str=sprintf('Threshold = %.2f',handles.Y(handles.IMAGE_NUMBER));
            set(handles.threshInt, 'String', str);
            set(handles.drawerSlider, 'Enable', 'off');
            set(handles.remove_single, 'Enable', 'off');
            set(handles.reset, 'Enable', 'off');
            set(handles.threshold_range, 'Enable', 'off');
            
            set(handles.smoo_text,'ForegroundColor',[.5 .5 .5]);
            set(handles.edit_smoo,'Enable','off','String','  ');
            
        otherwise
    
    end
    
    guidata(hObject, handles);
% -------------------------------------------------------------------------    



% --- Checkboxes ----------------------------------------------------------
function boost_Callback(hObject, eventdata, handles)
    if get(hObject,'Value')
        import(hObject,1);
    else
        import(hObject,0);
    end
% -------------------------------------------------------------------------    



% --- Editing -------------------------------------------------------------
function edit_smoo_Callback(hObject, eventdata, handles)

    smoo=round(abs(str2double(get(hObject,'String'))));
    
    if smoo > 1 && smoo < handles.NUM_ELEMENTS
        handles.SMOO = smoo;  
        
        axes(handles.histogramLoc);
        Y = GFilterA(handles.Y, handles.SMOO);
        X = handles.X;
        handles.DEFAULT_GAUSSIAN = Y;
        Y = Y-(handles.MEAN-handles.gaussianSavedSliderValue);
        handles.GAUSSIAN = Y;

        if ishandle(handles.CUST_LINE)
            if(handles.CUST_LINE ~= 0) 
                delete(handles.CUST_LINE);
            end
        end
        handles.CUST_LINE = line('Xdata', X, 'YData', Y, 'Color', 'r', 'EraseMode', 'normal','Marker','.'); 
        str=sprintf('Threshold = %.2f',handles.GAUSSIAN(handles.IMAGE_NUMBER));
        set(handles.threshInt, 'String', str);
        %set(handles.remove_single, 'Enable', 'off');
        %set(handles.drawerSlider, 'Enable', 'off');
        %set(handles.reset, 'Enable', 'off');
        %set(handles.threshold_range, 'Enable', 'off');
        
        
    else
        errordlg('Smoothing window makes sense only if it is a number larger that 1 and smaller than the number of frames');
    end
    set(hObject,'String',num2str(handles.SMOO));
    guidata(hObject, handles);

function edit_smoo_CreateFcn(hObject, eventdata, handles)

    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
% -------------------------------------------------------------------------    



% --- Manu ----------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)

function importMovie_Callback(hObject, eventdata, handles)

    set(handles.boost,'Value',0);
    import(hObject,0);
 
function save_as_Callback(hObject, eventdata, handles)

    dlg_title = 'Save Masked Images as a TIF-File';
    [FileName,PathName,FilterIndex] = uiputfile('*.tif',dlg_title,'MaskedImages');
    if FilterIndex
        if handles.LastSel==1
            thresh=handles.Y;
        elseif handles.LastSel==2
            thresh=handles.GAUSSIAN;
        elseif handles.LastSel==3
            thresh=handles.CUST_VALUES;
        end  
       
        axes(handles.progress1);
        cla;     
        bx=[0 1 1 0 0]/handles.NUM_ELEMENTS;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);     
        
        im=handles.IMAGECELL{1};
        
        %imt=im; imt(im<thresh(1))=0; imt(im>=thresh(1))=1;
        imt=im-thresh(1);
        imt(imt<0)=0;
        if max(imt(:))==min(imt(:))
            imt=ones(size(im));
        end
        
        [L,N]=bwlabel(imfill(logical(imt),'holes')); A=zeros(1,N);
        for a=1:N 
            A(a)=sum(sum(L==a));
        end
        [~,ma]=max(A);
        imt=(L==ma);
        masked=uint16(imt);
        try
            imwrite(masked,[PathName FileName],'Compression','none');            
        catch
            try
                imwrite(masked,[PathName FileName],'Compression','none','WriteMode','append');                    
            catch
                disp(['frame ' num2str(1) ' is NOT saved']);
            end  
        end
        for i=2:handles.NUM_ELEMENTS
            axes(handles.progress1);
            cla;     
            bx=[0 i i 0 0]/handles.NUM_ELEMENTS;
            by=[0 0 1 1 0];
            fill(bx,by,[.7 .7 .7]);
            set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
            pause(0.001);     
                
            im=handles.IMAGECELL{i};
            %imt=im; imt(im<thresh(i))=0; imt(im>=thresh(i))=1;
            imt=im-thresh(i);
            imt(imt<0)=0;
            if max(imt(:))==min(imt(:))
                imt=ones(size(im));
            end
            [L,N]=bwlabel(imfill(logical(imt),'holes')); A=zeros(1,N);           
            for a=1:N 
                A(a)=sum(sum(L==a));
            end
            [~,ma]=max(A);
            imt=(L==ma);
            masked=uint16(imt);
            status=0;
            for tr=1:10
                try
                    imwrite(masked,[PathName FileName],'Compression','none','WriteMode','append');                
                catch
                    pause(1);
                    continue;
                end
                status=1;
                break;
            end
            %disp([i tr])
            if ~status
                disp(['frame ' num2str(i) ' is NOT saved']);                                
            end
        end
        cla(handles.progress1,'reset');    
        set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
    end
%{
function single_file_Callback(hObject, eventdata, handles)

function sequence_Callback(hObject, eventdata, handles)

function seqn_bound_Callback(hObject, eventdata, handles)

    prompt = {'Name:','Start At:','Digits (1-8)'};
    dlg_title = 'Save File Sequence';
    num_lines = 1;
    def = {'BOUND','1','3'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        dig=round(abs(str2num(answer{3})));    
        if ~isempty(dig) && dig>0 && dig<9  
            str=['%0' answer{3} 'd'];
        else
            str='%03d';
        end
        beg=round(abs(str2num(answer{2}))); 
        if ~isempty(beg)
            fnum=sprintf(str,beg);
        else
            fnum=sprintf(str,1);
        end
        [FileName,PathName,FilterIndex] = uiputfile('*.txt',dlg_title,[answer{1} fnum]);

        if FilterIndex    
            for i=1:handles.NUM_ELEMENTS
                if mod(beg,1)==0
                    fnum=sprintf(str,beg+i-1);
                else
                    fnum=sprintf(str,i);
                end
                bnd=handles.BOUNDARIES{i};
                save([PathName answer{1} fnum '.txt'], 'bnd', '-ASCII');
            end
        end

    end

function seqn_thresh_Callback(hObject, eventdata, handles)
    prompt = {'Name:','Start At:','Digits (1-8)'};
    dlg_title = 'Save Image Sequence';
    num_lines = 1;
    def = {'THRESH','1','3'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        dig=round(abs(str2num(answer{3})));    
        if ~isempty(dig) && dig>0 && dig<9  
            str=['%0' answer{3} 'd'];
        else
            str='%03d';
        end
        beg=round(abs(str2num(answer{2}))); 
        if ~isempty(beg)
            fnum=sprintf(str,beg);
        else
            fnum=sprintf(str,1);
        end
        [FileName,PathName,FilterIndex] = uiputfile('*.tif',dlg_title,[answer{1} fnum]);

        if FilterIndex 
            if handles.LastSel==1
                thresh=handles.Y;
            elseif handles.LastSel==2
                thresh=handles.GAUSSIAN;
            elseif handles.LastSel==3
                thresh=handles.CUST_VALUES;
            end  
            for i=1:handles.NUM_ELEMENTS
                if mod(beg,1)==0
                    fnum=sprintf(str,beg+i-1);
                else
                    fnum=sprintf(str,i);
                end                
                im=handles.IMAGECELL{i};
                imt=im; imt(im<thresh(i))=0;
                thresholded=uint16(imt);
                imwrite(thresholded,[PathName answer{1} fnum '.tif'],'Compression','none');
            end
        end
    end

function seqn_mask_Callback(hObject, eventdata, handles)
    prompt = {'Name:','Start At:','Digits (1-8)'};
    dlg_title = 'Save Image Sequence';
    num_lines = 1;
    def = {'MASK','1','3'};
    answer = inputdlg(prompt,dlg_title,num_lines,def);

    if ~isempty(answer)
        dig=round(abs(str2num(answer{3})));    
        if ~isempty(dig) && dig>0 && dig<9  
            str=['%0' answer{3} 'd'];
        else
            str='%03d';
        end
        beg=round(abs(str2num(answer{2}))); 
        if ~isempty(beg)
            fnum=sprintf(str,beg);
        else
            fnum=sprintf(str,1);
        end
        [FileName,PathName,FilterIndex] = uiputfile('*.tif',dlg_title,[answer{1} fnum]);

        if FilterIndex    
            if handles.LastSel==1
                thresh=handles.Y;
            elseif handles.LastSel==2
                thresh=handles.GAUSSIAN;
            elseif handles.LastSel==3
                thresh=handles.CUST_VALUES;
            end  
            for i=1:handles.NUM_ELEMENTS
                if mod(beg,1)==0
                    fnum=sprintf(str,beg+i-1);
                else
                    fnum=sprintf(str,i);
                end                
                im=handles.IMAGECELL{i};
                imt=im; imt(im<thresh(i))=0; imt(im>=thresh(i))=1; 
                masked=uint16(imt);
                imwrite(masked,[PathName answer{1} fnum '.tif'],'Compression','none');
            end
        end
    end
    
function sing_bound_Callback(hObject, eventdata, handles)
    boundaries=handles.BOUNDARIES;
    uisave('boundaries');

function sing_thresh_Callback(hObject, eventdata, handles)
    dlg_title = 'Save Thresholded Images as a TIF-File';
    [FileName,PathName,FilterIndex] = uiputfile('*.tif',dlg_title,'ThreshImages');
    if FilterIndex
        if handles.LastSel==1
            thresh=handles.Y;
        elseif handles.LastSel==2
            thresh=handles.GAUSSIAN;
        elseif handles.LastSel==3
            thresh=handles.CUST_VALUES;
        end  
        im=handles.IMAGECELL{1};
        imt=im; imt(im<thresh(1))=0;
        thresholded=uint16(imt);
        imwrite(thresholded,[PathName FileName],'Compression','none');
        for i=2:handles.NUM_ELEMENTS
            im=handles.IMAGECELL{i};
            imt=im; imt(im<thresh(i))=0;
            thresholded=uint16(imt);
            imwrite(thresholded,[PathName FileName],'Compression','none','WriteMode','append');            
        end
    end
  

function sing_every_Callback(hObject, eventdata, handles)
% hObject    handle to sing_every (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    boundaries=handles.BOUNDARIES;
    N=handles.NUM_ELEMENTS;
    if handles.LastSel==1
        thresh=handles.Y;
    elseif handles.LastSel==2
        thresh=handles.GAUSSIAN;
    elseif handles.LastSel==3
        thresh=handles.CUST_VALUES;
    end    
    thresholded=cell(1,N);
    masked=cell(1,N);
    for i=1:N
        im=handles.IMAGECELL{i};
        imt=im; imt(im<thresh(i))=0;
        thresholded{i}=imt;
        imt(im>=thresh(i))=1;
        masked{i}=imt;
    end        
    uisave({'boundaries','N','thresh','thresholded','masked'});
  
function sing_mask_Callback(hObject, eventdata, handles)
    dlg_title = 'Save Masked Images as a TIF-File';
    [FileName,PathName,FilterIndex] = uiputfile('*.tif',dlg_title,'MaskedImages');
    if FilterIndex
        if handles.LastSel==1
            thresh=handles.Y;
        elseif handles.LastSel==2
            thresh=handles.GAUSSIAN;
        elseif handles.LastSel==3
            thresh=handles.CUST_VALUES;
        end  
        
        axes(handles.progress1);
        cla;     
        bx=[0 1 1 0 0]/handles.NUM_ELEMENTS;
        by=[0 0 1 1 0];
        fill(bx,by,[.7 .7 .7]);
        set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
        pause(0.001);     
        
        im=handles.IMAGECELL{1};
        imt=im; imt(im<thresh(1))=0; imt(im>=thresh(1))=1;
        [L,N]=bwlabel(imfill(imt,'holes')); A=zeros(1,N);
        for a=1:N 
            A(a)=sum(sum(L==a));
        end
        [~,ma]=max(A);
        imt=(L==ma);
        masked=uint16(imt);
        try
            imwrite(masked,[PathName FileName],'Compression','none');            
        catch
            try
                imwrite(masked,[PathName FileName],'Compression','none','WriteMode','append');                    
            catch
                disp(['frame ' num2str(1) ' is NOT saved']);
            end  
        end
        for i=2:handles.NUM_ELEMENTS
            axes(handles.progress1);
            cla;     
            bx=[0 i i 0 0]/handles.NUM_ELEMENTS;
            by=[0 0 1 1 0];
            fill(bx,by,[.7 .7 .7]);
            set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
            pause(0.001);     
                
            im=handles.IMAGECELL{i};
            imt=im; imt(im<thresh(i))=0; imt(im>=thresh(i))=1;
            [L,N]=bwlabel(imfill(imt,'holes')); A=zeros(1,N);
            for a=1:N 
                A(a)=sum(sum(L==a));
            end
            [~,ma]=max(A);
            imt=(L==ma);
            masked=uint16(imt);
            status=0;
            for tr=1:10
                try
                    imwrite(masked,[PathName FileName],'Compression','none','WriteMode','append');                
                catch
                    pause(1);
                    continue;
                end
                status=1;
                break;
            end
            %disp([i tr])
            if ~status
                disp(['frame ' num2str(i) ' is NOT saved']);                                
            end
        end
        cla(handles.progress1,'reset');    
        set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');
    end
%}    
function opt_Callback(hObject, eventdata, handles)

function cmap_Callback(hObject, eventdata, handles)

    str={'jet','gray','hot','cool','bone','hsv','spring','summer','autumn',...
         'winter','pink','copper','colorcube','prism'};

    [s,v] = listdlg('Name','Colormap',...
                    'PromptString','Select the map:',...
                    'SelectionMode','single',...
                    'ListSize',[160 200],...
                    'ListString',str);
    if v
        colormap(handles.imageLoc,str{s});
    end

function blinecol_Callback(hObject, eventdata, handles)

    str={'Yellow','Black','White','Red','Green','Blue','Cyan','Magenta'};
    col={'y','k','w','r','g','b','c','m'};
    [s,v] = listdlg('Name','Boundary Line',...
                    'PromptString','Select the color:',...
                    'SelectionMode','single',...
                    'ListSize',[160 120],...
                    'ListString',str);
    if v
        handles.BLCOL=col{s};

        axes(handles.imageLoc);
        xlim=get(handles.imageLoc,'XLim');
        ylim=get(handles.imageLoc,'YLim');
        cla;
        hold on;
        imagesc(handles.IMAGECELL{handles.IMAGE_NUMBER}, 'Parent', handles.imageLoc); 
        axis image;
        axis ij;
        set(handles.imageLoc,'XLim',xlim);
        set(handles.imageLoc,'YLim',ylim);
        axis off;
        bound=handles.BOUNDARIES{handles.IMAGE_NUMBER};
        by=bound(:,1); bx=bound(:,2);
        plot(bx,by,handles.BLCOL,'LineWidth',2);

    end
    guidata(hObject, handles);

function algo_Callback(hObject, eventdata, handles)
    str={'Histogram','Otsu (MATLAB)'};
    [s,v] = listdlg('Name','Algorithm',...
                    'PromptString','Select the algorithm:',...
                    'SelectionMode','single',...
                    'ListSize',[160 30],...
                    'ListString',str);
    if v       
        handles.Algorithm=s;        
    end
    guidata(hObject, handles);
    
    %delete(handles.lisdraw1);
    %delete(handles.lisdraw2);

    x = 1:handles.NUM_ELEMENTS;
    y = zeros(1,handles.NUM_ELEMENTS);
    b = cell(1,handles.NUM_ELEMENTS);

    
    for k = 1:handles.NUM_ELEMENTS
       axes(handles.progress1);
       cla;
       bx=[0 k k 0 0]/handles.NUM_ELEMENTS;
       by=[0 0 1 1 0];
       fill(bx,by,[.7 .7 .7]);
       set(handles.progress1,'Box','on','XLim',[0 1],'YLim',[0 1],'XTick',[],'YTick',[]);
       pause(0.001); 
       %set(handles.Frame,'String',['Frame # ' num2str(k)]); 
       [y(k), b{k}] = threshold(hObject, handles.IMAGECELL{k});             
    end
    cla(handles.progress1,'reset');   
    set(handles.progress1,'Box','on','XTick',[],'YTick',[],'Visible','off');

    handles.X = x;
    handles.Y = y;
    handles.BOUNDARIES = b;

    gaussianY = GFilterA(y, handles.SMOO); 
    handles.GAUSSIAN = gaussianY;

    handles.MEAN = mean(y);
    handles.CUST_VALUES = repmat(handles.MEAN, 1, handles.NUM_ELEMENTS);
    handles.CUST_X=[1 handles.NUM_ELEMENTS];
    handles.CUST_Y=[handles.MEAN handles.MEAN];
    
    axes(handles.imageLoc);
    cla;
    hold on;
    imagesc(handles.IMAGECELL{1});
    bound=b{1};
    by=bound(:,1); bx=bound(:,2);
    plot(bx,by,handles.BLCOL,'LineWidth',2);
    axis image; 
    axis ij;
    axis off;

    if handles.NUM_ELEMENTS > 1
        axes(handles.histogramLoc);
        cla;
        hold on;
        plot(x, y, '.-'); 
        range=max(y)-min(y);
        if range>0
            set(handles.histogramLoc, 'YLim', [floor(min(y)-range/3) ceil(max(y)+range/3)]);
        else                
            set(handles.histogramLoc, 'YLim', [floor(max(y)-1) ceil(max(y)+1)]);
        end       

        if ~ishandle(handles.REDLINE) || handles.REDLINE==0            
            handles.REDLINE = line([1, 1], get(handles.histogramLoc, 'YLim'), 'Color', 'r', 'Parent', ...
                handles.histogramLoc, 'EraseMode', 'normal'); 
        end
        
        if range>0
            set(handles.drawerSlider, 'Min', floor(min(y)-range/3), 'Max', ceil(max(y)+range/3));
        else                
            set(handles.drawerSlider, 'Min', floor(max(y)-1), 'Max', ceil(max(y)+1));
        end
    else
        set(handles.drawerSlider, 'Min', max(floor(0.9*max(y)),0), 'Max', max(ceil(1.1*max(y)),1));
        
    end

    guidata(hObject, handles);
    set(handles.drawerSlider, 'Value',handles.MEAN);
    set(handles.movieSlider, 'Value',1);
   
    set(handles.suggested, 'Value', 1);
    set(handles.remove_single, 'Enable', 'off');
    set(handles.reset, 'Enable', 'off');

    str=sprintf('Threshold = %.2f',y(1));
    set(handles.threshInt, 'String', str);
    set(handles.Frame,'String','Frame # 1');
                 
% -------------------------------------------------------------------------    

function set_initial_min_thresh(val)
    handles=guidata(gcf);
    set(handles.min_thresh_range, 'String', num2str(val));
    handles.MIN_THRESH_RANGE = val;
    guidata(gcf, handles);

function set_initial_max_thresh(val)
    handles=guidata(gcf);
    set(handles.max_thresh_range, 'String', num2str(val));
    handles.MAX_THRESH_RANGE = val;
    disp(val);
    guidata(gcf, handles);
        
        
function update_thresh_range()
% --- Allows user to customize the upper and lower bounds of the threshold
% range. The range of the intensity plot will change as a result and any
% thresholding of a single image will be done within these bounds.

    handles = guidata(gcf);

    answer = { get(handles.max_thresh_range, 'String'), get(handles.min_thresh_range, 'String') };
    
    if ~isempty(answer)
        % validate
        up=ceil(abs(str2num(answer{1})));
        lw=ceil(abs(str2num(answer{2})));
        
        if length(up) == 1 && length(lw) == 1 && up > lw
 
            v=get(handles.drawerSlider,'Value');
            handles.ON_CURRENT_FRAME = 0;

            if handles.NUM_ELEMENTS>1

                ma=max(handles.CUST_Y);
                mi=min(handles.CUST_Y);
                
                if ~isempty(up) && ~isempty(lw) 
                    if up>=ma && lw<=mi 
                                                                       
                        set(handles.histogramLoc,'YLim',[lw up]);
                        handles.MAX_THRESH_RANGE = up;
                        handles.MIN_THRESH_RANGE = lw;
                        if v>up
                            guidata(gcf, handles);
                            set(handles.drawerSlider,'Value',up);
                            handles = guidata(gcf);
                            str=sprintf('Threshold = %.2f',up);
                            set(handles.threshInt, 'String', str);
                        elseif v<lw
                            guidata(gcf, handles);
                            set(handles.drawerSlider,'Value',lw);
                            handles = guidata(gcf);
                            str=sprintf('Threshold = %.2f',lw);
                            set(handles.threshInt, 'String', str);
                        end
                        set(handles.drawerSlider,'Min',lw,'Max',up);
                        set(handles.min_thresh_range, 'String', num2str(lw));
                        set(handles.max_thresh_range, 'String', num2str(up));

                        if ishandle(handles.REDLINE) && handles.REDLINE ~= 0
                            delete(handles.REDLINE);
                            handles.REDLINE = line([handles.IMAGE_NUMBER, handles.IMAGE_NUMBER], get(handles.histogramLoc, 'YLim'), 'Color', 'r', 'Parent', ...
                            handles.histogramLoc, 'EraseMode', 'normal'); 
                        end
                    else
                        set(handles.max_thresh_range, 'String', num2str(handles.MAX_THRESH_RANGE));
                        set(handles.min_thresh_range, 'String', num2str(handles.MIN_THRESH_RANGE));
                    end            
                end
            elseif handles.NUM_ELEMENTS==1
                set(handles.drawerSlider,'Enable','off');
                if ~isempty(up) && ~isempty(lw)
                    handles.MAX_THRESH_RANGE = up;
                    handles.MIN_THRESH_RANGE = lw;
                    if v>up
                        set(handles.drawerSlider,'Value',up);
                        str=sprintf('Threshold = %.2f',up);
                        set(handles.threshInt, 'String', str);
                    elseif v<lw
                        set(handles.drawerSlider,'Value',lw);
                        str=sprintf('Threshold = %.2f',lw);
                        set(handles.threshInt, 'String', str);
                    end
                    set(handles.drawerSlider,'Min',lw,'Max',up);
                    set(handles.min_thresh_range, 'String', num2str(lw));
                    set(handles.max_thresh_range, 'String', num2str(up));
                else
                    set(handles.max_thresh_range, 'String', num2str(handles.MAX_THRESH_RANGE));
                    set(handles.min_thresh_range, 'String', num2str(handles.MIN_THRESH_RANGE));
                end            
            end
        else
            set(handles.max_thresh_range, 'String', num2str(handles.MAX_THRESH_RANGE));
            set(handles.min_thresh_range, 'String', num2str(handles.MIN_THRESH_RANGE));
        end
    end
    
    

    guidata(gcf, handles);


function min_thresh_range_Callback(hObject, eventdata, handles)
update_thresh_range();

% --- Executes during object creation, after setting all properties.
function min_thresh_range_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function max_thresh_range_Callback(hObject, eventdata, handles)
update_thresh_range();

% --- Executes during object creation, after setting all properties.
function max_thresh_range_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
