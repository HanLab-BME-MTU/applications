function editPropertiesGUI_OpeningFcn(hObject, eventdata, handles, varargin)
%--executes just before editPropertiesGUI is made visible

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to editPropertiesGUI (see VARARGIN)

% Choose default command line output for editPropertiesGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


amgH = findobj('Tag','AMG');
amgHandles = guidata(amgH);
handles.amgHandles = amgHandles;

%check if any job loaded
dirList = get(amgHandles.dirListBox,'String');
if ~iscell(dirList)
    close(hObject);
    return
    %no cell, no job, no execution
end

%check if default biodata-dir exists
mainDir = getenv('BIODATA');
if isempty(mainDir)
    mainDir = getenv('HOME');
end

activeJobNum = get(amgHandles.dirListBox,'Value');
activeJobNameList = get(amgHandles.dirListBox,'String');
activeJobName = activeJobNameList(activeJobNum);
myJob = amgHandles.job(activeJobNum);
handles.activeJobNum = activeJobNum;

%set gui title
set(handles.EPGUI,'Name',['Edit properties for: ',char(activeJobName)]);

%get the checkbox handles and store them in a vector
handles.checkH = [handles.edit_check_filter;...
        handles.edit_check_detect;...
        handles.edit_check_link;...
        handles.edit_check_labelgui;...
        handles.edit_check_track;...
        handles.edit_check_labelgui2;...
        handles.edit_check_analysis1;...
        handles.edit_check_analysis2];

handles.cellCycleH = [handles.edit_check_CC_G1;...
        handles.edit_check_CC_Start;...
        handles.edit_check_CC_SPhase;...
        handles.edit_check_CC_G2;...
        handles.edit_check_CC_Meta;...
        handles.edit_check_CC_Ana;...
        handles.edit_check_CC_Telo];

handles.strainH = [handles.edit_check_strain_WT;...
        handles.edit_check_strain_dam1_1;...
        handles.edit_check_strain_dam1_11;...
        handles.edit_check_strain_ipl1_321;...
        handles.edit_check_strain_ipl1_1;...
        handles.edit_check_strain_ndc10_1;...
        handles.edit_check_strain_ndc80_1;...
        handles.edit_check_strain_bik1;...
        handles.edit_check_strain_bim1];

handles.drugH = [handles.edit_check_drugs_nocodazole;...
        handles.edit_check_drugs_benomyl;...
        handles.edit_check_drugs_hydroxyurea;...
        handles.edit_check_drugs_alphaFactor];

handles.imgH = [handles.edit_check_split;...
        handles.edit_timeptsInMem_txt];

handles.subsetH = [handles.edit_check_analyzeSelected;...
        handles.edit_startTime_txt;...
        handles.edit_endTime_txt];

handles.detectH = [handles.edit_ftest_prob_txt;...
        handles.edit_maxslope_txt];

handles.linkH = [handles.edit_IDopt_checkIntensity_PD;...
        handles.edit_IDopt_verbose_PD;...
        handles.edit_IDopt_weight_txt];

handles.backgroundH = [handles.edit_check_correctBG;...
        handles.edit_RB_standardBG;...
        handles.edit_standardMovie1_PD;...
        handles.edit_standardMovie2_PD;...
        handles.edit_RB_firstAndLast;...
        handles.edit_txt_useFirst;...
        handles.edit_txt_useLast;...
        handles.bgText1;...
        handles.bgText2];
        

%get list of possible temperatures
tempList = get(handles.edit_temperature_PD,'String');

%get as much data as possible:
%1. read r3dheader for microscopy data
%2. look for dataProperties in the job file or if it doesn't exist: look for a saveFile

%read r3d extended header
cd(mainDir);
cd(myJob.projProperties.dataPath);
if strcmp(myJob.fileExtension(end),'c')
    load('r3dMovieHeader.mat');
    if ~exist('r3dMovieHeader','var')
        error('can''t find movie header')
    end
    header = r3dMovieHeader;
    movieName = [myJob.projName,'.r3c'];
else
    movieName = [myJob.projName,'.r3d'];
    header = readr3dheader(movieName);
    header.cropInfo = [];
    header.correctInfo = [];
end

%store header
handles.header = header;

%check wavelength
if length(header.wvl) > 1
    close(hObject)
    errordlg('Sorry, multiple wavelengths not supported yet')
    error('Sorry, multiple wavelengths not supported yet')
    return
end

%check if there is already a filtered movie
filteredMovieName = chooseFile(['filtered_movie_',myJob.projName],[],'new');
if ~isempty(filteredMovieName)
    if myJob.projProperties.status==0
        myJob.projProperties.status = 1;
    end
end


%-----------------------set data----------------------------------------


%movie data
pixelX = sprintf('%0.4f',header.pixelX);
pixelZ = sprintf('%0.4f',header.pixelZ);

set(handles.edit_pixel_xy_txt,'String',pixelX);
set(handles.edit_pixel_z_txt,'String',pixelZ);
set(handles.edit_lensID_txt,'String',header.lensID);
if header.lensID==12003 %do like this until we can look it up somewhere
    NA = 1.4;
    set(handles.edit_NA_txt,'String',1.4); 
elseif header.lensID==10002 %do like this until we can look it up somewhere
    NA = 1.4;
    set(handles.edit_NA_txt,'String',1.4); 
    
else
    NA = [];
    set(handles.edit_NA_txt,'String','','Style','edit');
end
set(handles.edit_wavelength_txt,'String',header.wvl(1));
set(handles.edit_exposureTime_txt,'String',header.expTime);
set(handles.edit_movieLength_txt,'String',header.numTimepoints);
set(handles.edit_movieSize_x_txt,'String',header.numRows);
set(handles.edit_movieSize_y_txt,'String',header.numCols);
set(handles.edit_movieSize_z_txt,'String',header.numZSlices);
set(handles.edit_NDfilter_txt,'String',header.ndFilter);

%get frameTime and meanTimeLapse
frameTime = header.Time;
frameTime = reshape(frameTime,header.numZSlices,header.numTimepoints);
deltaFrameTime = frameTime(:,2:end)-frameTime(:,1:end-1);
if ~isempty(deltaFrameTime) %if there's only one frame, mean([]) will produce a warning
    meanTimeLapse = mean(deltaFrameTime(:));
else
    meanTimeLapse = NaN;
end
timeLapse = sprintf('%0.3f',meanTimeLapse);
set(handles.edit_timeLapse_txt,'String',timeLapse);

%filter parameters (can't be edited)
%start assuming that the psf matches the theoretical psf
handles.sigmaCorrection(1) = 1.12;%1.2 before pixCorr.
handles.sigmaCorrection(2) = 1.0;%1.3 before NAcorr

set(handles.edit_sigmaCorrX_txt,'String',handles.sigmaCorrection(1));
set(handles.edit_sigmaCorrY_txt,'String',handles.sigmaCorrection(1));
set(handles.edit_sigmaCorrZ_txt,'String',handles.sigmaCorrection(2));

if ~isempty(NA)
    [FT_XY, FT_Z] = calcFilterParms(...
        header.wvl(1),NA,1.51,'gauss',handles.sigmaCorrection, [header.pixelX, header.pixelZ]);

    patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
    handles.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];
    set(handles.edit_filterX_txt,'String',sprintf('%0.3f',FT_XY));
    set(handles.edit_filterY_txt,'String',sprintf('%0.3f',FT_XY));
    set(handles.edit_filterZ_txt,'String',sprintf('%0.3f',FT_Z));
    set(handles.edit_patchX_txt,'String',patchXYZ(1));
    set(handles.edit_patchY_txt,'String',patchXYZ(2));
    set(handles.edit_patchZ_txt,'String',patchXYZ(3));
end

%store data
handles.frameTime = frameTime';
handles.meanTimeLapse = meanTimeLapse;
handles.pixelsizeXYZ = [header.pixelX,header.pixelZ];


%find (if exists) project-data
projData = chooseFile([myJob.projName,'-data'],[],'GUI','log');


if ~isempty(myJob.dataProperties)&isfield(myJob.dataProperties,'IDopt')
    
    %if several projData and user cancelled: take job-projData
    if isempty(projData)
        handles.projData = myJob.projData;
        createNew = 1; %if the user does not want to load projData, he/she probably wants to redo everything
    else
        handles.projData = projData;
        createNew = 0;
    end
    
    %there is already data for this project
    if isfield(myJob.dataProperties,'cellCycle')
        set(handles.cellCycleH(find(bsum2lvec(myJob.dataProperties.cellCycle))),'Value',1);
        set(handles.strainH(find(bsum2lvec(myJob.dataProperties.strains))),'Value',1);
        set(handles.drugH(find(bsum2lvec(myJob.dataProperties.drugs))),'Value',1);
        set(handles.edit_temperature_PD,'Value',strmatch(myJob.dataProperties.temperature,tempList));
    end
    
    %     if isfield(myJob.dataProperties,'analyzeSubset')
    %         if myJob.dataProperties.analyzeSubset(1)~=1|myJob.dataProperties.analyzeSubset(2)~=header.numTimepoints
    %             set(handles.subsetH(1),'Value',1);
    %         end
    
    if myJob.projProperties.status>0
        setGreen = bsum2lvec(myJob.projProperties.status);
        set(handles.checkH(find(setGreen)),'BackgroundColor',[0,1,0]);
        set(handles.checkH,'Value',0);
    end
    
    %store status
    handles.status = myJob.projProperties.status;
    
    %mark jobs2run
    if isfield(myJob,'jobs2run')
        jobs2run = bsum2lvec(myJob.jobs2run);
        set(handles.checkH(find(jobs2run)),'Value',1);
    end
    
    %set eraseAllPrev
    if isfield(myJob,'eraseAllPrev')
        if ~isempty(myJob.eraseAllPrev)
            set(handles.edit_check_eraseAllPrevData,'Value',myJob.eraseAllPrev);
        end
    end
    
    %make sure createNew is not empty
    if isfield(myJob,'createNew')
        if isempty(myJob.createNew) %if it's not empty, user already has decided on something
            myJob.createNew = createNew;
        end
    else
        myJob.createNew = createNew;
    end
    
    %set autospfinder properties
    set(handles.edit_maxslope_txt,'String',myJob.dataProperties.CH_MAXSLOPE); %   1 pt/0.5 grayvalue (max. slope)
    set(handles.edit_ftest_prob_txt,'String',myJob.dataProperties.F_TEST_PROB);    % min probability that there are two spots in a distribution
    
    %set spotID properties
    if isfield(myJob.dataProperties,'IDopt') %prevent conflict with older versions
        set(handles.edit_IDopt_weight_txt,'String',myJob.dataProperties.IDopt.weight);
        set(handles.edit_IDopt_checkIntensity_PD,'Value',2-myJob.dataProperties.IDopt.checkIntensity);
        set(handles.edit_IDopt_verbose_PD,'Value',myJob.dataProperties.IDopt.verbose+1);
    end
    %
    % more 2 come
    %
else
    handles.projData = projData;
    if ~isempty(projData)
        %there is a property file
        load([projData],'dataProperties');
        myJob.dataProperties = dataProperties;
        load([projData],'projProperties');
        myJob.projProperties = projProperties;
        
        %if this file exists, the respective fields exist, too
        if isfield(myJob.dataProperties,'cellCycle')
            set(handles.cellCycleH(find(bsum2lvec(myJob.dataProperties.cellCycle))),'Value',1);
            set(handles.strainH(find(bsum2lvec(myJob.dataProperties.strains))),'Value',1);
            set(handles.drugH(find(bsum2lvec(myJob.dataProperties.drugs))),'Value',1);
            set(handles.edit_temperature_PD,'Value',strmatch(myJob.dataProperties.temperature,tempList));
        end
        
        
        
        %mark jobs already done
        if projProperties.status==0&myJob.projProperties.status==1
            projProperties.status = 1;
        end
        if projProperties.status>0
            setGreen = bsum2lvec(projProperties.status);
            set(handles.checkH(find(setGreen)),'BackgroundColor',[0,1,0]);
            set(handles.checkH,'Value',0);
        end
        
        %mark jobs2run
        if isfield(myJob,'jobs2run')
            jobs2run = bsum2lvec(myJob.jobs2run);
            set(handles.checkH(find(jobs2run)),'Value',1);
        end
        
        %set eraseAllPrev
        if isfield(myJob,'eraseAllPrev')
            if isempty(myJob.eraseAllPrev)
                myJob.eraseAllPrev = 0;
            end
            set(handles.edit_check_eraseAllPrevData,'Value',myJob.eraseAllPrev);
        end
        
        handles.status = projProperties.status;
        myJob.createNew = 0;
        
        %set autospfinder properties
        set(handles.edit_maxslope_txt,'String',dataProperties.CH_MAXSLOPE); %   1 pt/0.5 grayvalue (max. slope)
        set(handles.edit_ftest_prob_txt,'String',dataProperties.F_TEST_PROB);    % min probability that there are two spots in a distribution
        
        %set spotID properties
        if isfield(dataProperties,'IDopt')
            set(handles.edit_IDopt_weight_txt,'String',dataProperties.IDopt.weight);
            set(handles.edit_IDopt_checkIntensity_PD,'Value',2-dataProperties.IDopt.checkIntensity);
            set(handles.edit_IDopt_verbose_PD,'Value',dataProperties.IDopt.verbose+1);
        end
        %
        %set other properties, too
        %
        
    else 
        myJob.createNew = 1; %if no dataProperties-file: create new anyway
        handles.status = myJob.projProperties.status;
        %mark jobs already done
        if handles.status==1
            set(handles.checkH(1),'BackgroundColor',[0,1,0]);
            set(handles.checkH,'Value',0);
        end
        
        %if a movie has been "edited" before, there should be new defaults: set
        if isfield(amgHandles,'epguiDefaults')
            set(handles.cellCycleH,{'Value'},amgHandles.epguiDefaults.movieProperties.cellCycle);
            set(handles.strainH,{'Value'},amgHandles.epguiDefaults.movieProperties.strain);
            set(handles.drugH,{'Value'},amgHandles.epguiDefaults.movieProperties.drugs);
            set(handles.edit_temperature_PD,'Value',amgHandles.epguiDefaults.movieProperties.tempValue);
        end
    end
end

%init crop-field and set crop status
myJob.dataProperties.crop = handles.header.cropInfo; %empty, if r3dreadheader was called (r3d-movie)

if ~isempty(myJob.dataProperties.crop)
    set(handles.edit_movieSize_x_txt,'String',myJob.dataProperties.crop(2,1),'ForegroundColor','r');
    set(handles.edit_movieSize_y_txt,'String',myJob.dataProperties.crop(2,2),'ForegroundColor','r');
end

%prepare BG correction
handles.currentBGState = [];
handles.oldBGState = [];
handles.correctBackground1 = {};
handles.correctBackground2 = {};

if ~isfield(handles.header,'correctInfo')
    handles.header.correctInfo = [];
end

%fill pdList anyway (this does no harm)
if isfield(amgHandles,'epguiDefaults')
    handles.correctBackground1= amgHandles.epguiDefaults.background.hList1;
    handles.correctBackground2 = amgHandles.epguiDefaults.background.hList2;
    set(handles.backgroundH(3),'String',amgHandles.epguiDefaults.background.pdList1);
    set(handles.backgroundH(4),'String',amgHandles.epguiDefaults.background.pdList2);
    set(handles.backgroundH(6),'String',amgHandles.epguiDefaults.background.framesS);
    set(handles.backgroundH(7),'String',amgHandles.epguiDefaults.background.framesL);
end

%set backgroundCorrection-info (ci)
if ~isempty(handles.header.correctInfo)
    ci = handles.header.correctInfo;
elseif ~isempty(myJob.correctBackground)
    ci = myJob.correctBackground;
    set(handles.backgroundH(1),'Value',1)
    set(handles.backgroundH([2,5]),'Enable','on')
    if ~isempty(ci.correctFiles)
        set(handles.backgroundH([3:4]),'Enable','on');
    else
        set(handles.backgroundH([6:end]),'Enable','on');
    end
else
    ci = [];
end
if ~isempty(ci)
    if ~isempty(ci.correctFrames)
        set(handles.backgroundH(6),'String',ci.correctFrames(1));
        set(handles.backgroundH(7),'String',ci.correctFrames(2));
        set(handles.backgroundH(5),'Value',1);
        set(handles.backgroundH(2),'Value',0);
    elseif ~isempty(ci.correctFiles)
        pdString1 = get(handles.backgroundH(3),'String');
        cFileName = ci.correctFiles{1,2};
        pdString1 = [pdString1;{cFileName(1:end-4)}];
        v1 = length(pdString1);
        handles.correctBackground1 = [ci.correctFiles(1,:)];
        pdString2 = get(handles.backgroundH(4),'String');
        if size(ci.correctFiles,1)==2;
            cFileName = ci.correctFiles{2,2};
            pdString2 = [pdString2;{cFileName(1:end-4)}];
            v2 = length(pdString2);
            handles.correctBackground2 = [ci.correctFiles(2,:)];
        else
            v2 = 1;
        end
        set(handles.backgroundH(3:4),{'String'},{pdString1;pdString2},{'Value'},{v1;v2});
        set(handles.backgroundH(5),'Value',0);
        set(handles.backgroundH(2),'Value',1);
    end %~isempty(ci.correctFrames)
    
    ci.header = handles.header;
    handles.currentBGState = ci;
    handles.oldBGState = ci;
end %if ~isempty(ci)

%do not allow editing for jobs not being run
set(handles.edit_check_createNew,'Value',myJob.createNew);
if ~myJob.createNew
    if myJob.projProperties.status>0
        statVec = bsum2bvec(myJob.projProperties.status);
        
        if any(statVec==2)
            set(handles.detectH,'Enable','off');
        end
        
        if any(statVec==4)
            set(handles.linkH,'Enable','off');
        end
        
        %if there is a filtered movie: do not allow change of background, but set values
        %(if exist)
        if any(statVec == 1)
            set(handles.backgroundH,'Enable','off');
            
        end %if any(statVec == 1)
    end %if myJob.projProperties.status>0
end %~myJob.createNew

%if it is a corrected movie: do not allow another correction
if ~isempty(handles.header.correctInfo)
    set(handles.backgroundH,'Enable','off');
end


%check if additional PB's are allowed
%do not show both, as in updating jobs, apply2all is bugged
if amgHandles.acceptAll==0
    set(handles.edit_accept_all_PB,'Visible','on');
else
    if amgHandles.apply2all==0
        set(handles.edit_apply2all_PB,'Visible','on');
    end
end

%prepare preview
handles.previewData.filteredMovieName = filteredMovieName;
handles.previewData.movieName = movieName;
handles.previewData.movieLength = header.numTimepoints;


guidata(hObject,handles);
