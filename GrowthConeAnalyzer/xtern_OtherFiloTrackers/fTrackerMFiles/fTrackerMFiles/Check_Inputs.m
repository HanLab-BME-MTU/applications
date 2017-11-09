function Err = Check_Inputs(hObject, eventdata, handles)

% Author: Antoine Godin
% godin.antoine@sympatico.ca

    Err = 0;
    imagesFolder = get(handles.imagesFolder,'String');
    try 
        here = cd;
        cd(imagesFolder)
        cd(here)
        if length(dir([imagesFolder,'*.tif'])) == 0
            errordlg('No Images in given folder','Wrong Folder!');
            set(handles.imagesFolder,'Foregroundcolor',[1,0,0])
            Err = 1;
            return
        else
            set(handles.imagesFolder,'Foregroundcolor',[0,0,0])
        end
    catch
        errordlg('Problems with given folder','Wrong Folder!');
        set(handles.imagesFolder,'Foregroundcolor',[1,0,0])
        Err = 1;
        return
    end
            
    if str2num(get(handles.ThresholdFactor,'String')) < 0
        errordlg('Problems with Threshold Factor','Problems with Images Parameters!');
        set(handles.ThresholdFactor,'Foregroundcolor',[1,0,0])
        Err = 1;
        return
    else
        set(handles.ThresholdFactor,'Foregroundcolor',[0,0,0])
    end
                   
    if str2num(get(handles.PixelsDilate,'String')) < 0
        errordlg('Problems with Pixels Dilate','Problems with Images Parameters!');
        set(handles.PixelsDilate,'Foregroundcolor',[1,0,0])
        Err = 1;
        return
    else
        set(handles.PixelsDilate,'Foregroundcolor',[0,0,0])
    end

    if str2num(get(handles.DiskSize,'String')) < 0
        errordlg('Problems with Disk Size, it has to be an integer greater than 0','Problems with Images Parameters!');
        set(handles.DiskSize,'Foregroundcolor',[1,0,0])
        Err = 1;
        return
    else
        set(handles.DiskSize,'Foregroundcolor',[0,0,0])
    end

    if str2num(get(handles.CloseSize,'String')) < 0 
        errordlg('Problems with Close Size, it has to be a positive integer','Problems with Images Parameters!');
        set(handles.CloseSize,'Foregroundcolor',[1,0,0])
        Err = 1;
        return
    else
        set(handles.CloseSize,'Foregroundcolor',[0,0,0])
    end

    if str2num(get(handles.frameToTest,'String')) < get(handles.FrameToTestSlider,'Min') | str2num(get(handles.frameToTest,'String')) > get(handles.FrameToTestSlider,'Max')
        errordlg('Problems with Frame to Test','Problems with Image to Test!');
        set(handles.frameToTest,'Foregroundcolor',[1,0,0])
        Err = 1;
        return
    else
        set(handles.frameToTest,'Foregroundcolor',[0,0,0])
    end