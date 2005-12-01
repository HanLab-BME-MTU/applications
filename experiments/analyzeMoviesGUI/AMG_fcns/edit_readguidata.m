function edit_readguidata(handles)
%reads out data from editPropertiesGUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amgHandles = handles.amgHandles;
activeJobNum = handles.activeJobNum;

%read microscope data
%for sigma-calculation
PIXELSIZE_XY = handles.pixelsizeXYZ(1);
PIXELSIZE_Z = handles.pixelsizeXYZ(2);

NA = str2double(get(handles.edit_NA_txt,'String'));
if isempty(NA)
    errordlg('Specify NA first!');
    return
end
WVL = str2double(get(handles.edit_wavelength_txt,'String'));
%for storage
amgHandles.job(activeJobNum).dataProperties.PIXELSIZE_XY = PIXELSIZE_XY;
amgHandles.job(activeJobNum).dataProperties.PIXELSIZE_Z = PIXELSIZE_Z;
amgHandles.job(activeJobNum).dataProperties.LENSID = str2double(get(handles.edit_lensID_txt,'String'));
amgHandles.job(activeJobNum).dataProperties.NA = NA; %set until we know from lensID;
amgHandles.job(activeJobNum).dataProperties.WVL = WVL;
amgHandles.job(activeJobNum).dataProperties.timeLapse = handles.meanTimeLapse;
amgHandles.job(activeJobNum).dataProperties.frameTime = handles.frameTime;
amgHandles.job(activeJobNum).dataProperties.expTime = str2double(get(handles.edit_exposureTime_txt,'String'));
amgHandles.job(activeJobNum).dataProperties.movieSize = [str2double(get(handles.edit_movieSize_x_txt,'String')),...
        str2double(get(handles.edit_movieSize_y_txt,'String')),...
        str2double(get(handles.edit_movieSize_z_txt,'String')),...
        str2double(get(handles.edit_movieLength_txt,'String'))];
amgHandles.job(activeJobNum).dataProperties.NDfilter = str2double(get(handles.edit_NDfilter_txt,'String'));

%read movie properties data
ccCell = get(handles.cellCycleH,'Value');
amgHandles.job(activeJobNum).dataProperties.cellCycle = lvec2bsum(cat(1,ccCell{1:end}));
strainCell = get(handles.strainH,'Value');
amgHandles.job(activeJobNum).dataProperties.strains = lvec2bsum(cat(1,strainCell{1:end}));
drugCell = get(handles.drugH,'Value');
amgHandles.job(activeJobNum).dataProperties.drugs = lvec2bsum(cat(1,drugCell{1:end}));
tempList = get(handles.edit_temperature_PD,'String');
amgHandles.job(activeJobNum).dataProperties.temperature = tempList(get(handles.edit_temperature_PD,'Value'));

%write cropInfo (isempty if no cropping)
amgHandles.job(activeJobNum).dataProperties.crop = handles.header.cropInfo;

% maxSize. Multiply by 2e20 to get bytes 
maxSize = str2double(get(handles.edit_maxSize_txt,'String'))*2^20;
amgHandles.job(activeJobNum).dataProperties.maxSize = maxSize;

%set detector properties
MAXSPOTS = str2double(get(handles.edit_maxNumSpots_txt,'String'));
amgHandles.job(activeJobNum).dataProperties.MAXSPOTS = MAXSPOTS; 
F_TEST_PROB = str2double(get(handles.edit_ftest_prob_txt,'String'));
amgHandles.job(activeJobNum).dataProperties.F_TEST_PROB = F_TEST_PROB;    % min probability that there are two spots in a distribution

%set spotID properties
opt.weight = str2double(get(handles.edit_IDopt_weight_txt,'String'));
opt.checkIntensity = 2-get(handles.edit_IDopt_checkIntensity_PD,'Value');
opt.verbose = get(handles.edit_IDopt_verbose_PD,'Value')-1;
opt.save = 1; %you can not change that in the GUI
amgHandles.job(activeJobNum).dataProperties.IDopt = opt;


%write jobs to do
checkStatus = cell2mat(get(handles.checkH,'Value'));
jobs = lvec2bsum(checkStatus);
amgHandles.job(activeJobNum).jobs2run = jobs;

%write status of analysis
amgHandles.job(activeJobNum).projProperties.status = handles.status;

%write correction status
%get currentBGState
if get(handles.edit_check_correctBG,'Value')
    currentBGState = handles.currentBGState;
    %check that BGState is ok
    if isempty(currentBGState.correctFiles)
        if isempty(currentBGState.correctFrames)
            %there is no valid BGInfo at all
            errordlg('if you want to do background subtraction, please specify correct settings')
            return
        else %we'r correcting via frames, ergo dataProperties.movieSize and frameTime have to be adjusted
            %get frameTime and meanTimeLapse
            frameTime = handles.header.Time; %header will be changed when correctedMovie is calculated
            frameTime = reshape(frameTime,handles.header.numZSlices,handles.header.numTimepoints);
            %delete surplus frames
            frameTime = frameTime(:,currentBGState.correctFrames(1)+1:end-currentBGState.correctFrames(2));
            deltaFrameTime = diff(frameTime,1,2);
            if ~isempty(deltaFrameTime) %if there's only one frame, mean([]) will produce a warning
                meanTimeLapse = mean(deltaFrameTime(:));
            else
                meanTimeLapse = NaN;
            end
            amgHandles.job(activeJobNum).dataProperties.timeLapse = meanTimeLapse;
            amgHandles.job(activeJobNum).dataProperties.frameTime = frameTime';
        end
    end
    amgHandles.job(activeJobNum).correctBackground = currentBGState;
else
    amgHandles.job(activeJobNum).correctBackground = [];
end

%correct for movies starting at later frames only - they should start with t=0
amgHandles.job(activeJobNum).dataProperties.frameTime = ...
    amgHandles.job(activeJobNum).dataProperties.frameTime - min(amgHandles.job(activeJobNum).dataProperties.frameTime(:));

%write all the not-yet and never-to-be implemented variables into dataProperties
%never-to-be
amgHandles.job(activeJobNum).dataProperties.PATCHSIZE = 7;
amgHandles.job(activeJobNum).dataProperties.CH_MAXNUMINTERV =  1000;
amgHandles.job(activeJobNum).dataProperties.OVERLPSIZE =  [15 15 15];

%filter parameters
amgHandles.job(activeJobNum).dataProperties.sigmaCorrection = handles.sigmaCorrection;
amgHandles.job(activeJobNum).dataProperties.FT_SIGMA     = handles.FILTERPRM(1:3);
amgHandles.job(activeJobNum).dataProperties.FILTERPRM = handles.FILTERPRM;

%not-yet
amgHandles.job(activeJobNum).dataProperties.split.set = 0;
amgHandles.job(activeJobNum).dataProperties.split.TIMEPTS_IN_MEM =  10;

amgHandles.job(activeJobNum).dataProperties.T_TEST_PROB =  0.05;    % min prob. that distance is zero

amgHandles.job(activeJobNum).projData = handles.projData; %name of data-file
amgHandles.job(activeJobNum).createNew = get(handles.edit_check_createNew,'Value'); %if 0, data is appended to existing project

amgHandles.job(activeJobNum).eraseAllPrev = get(handles.edit_check_eraseAllPrevData,'Value'); %erases all older data files

%write name of file/project
amgHandles.job(activeJobNum).projName = handles.previewData.movieName(1:end-4);
%amgHandles.job(activeJobNum).fileExtension = handles.previewData.movieName(end-3:end);
amgHandles.job(activeJobNum).dataProperties.name = amgHandles.job(activeJobNum).projName;

%update amg-list-entry
dirList = get(amgHandles.dirListBox,'String');
dirList{activeJobNum} = amgHandles.job(activeJobNum).projName;
%store new dirList
set(amgHandles.dirListBox,'String',dirList);

%write help
amgHandles.job(activeJobNum).dataProperties.help.status = [...
        {'project status'};...
        {'2: spot detector'};...
        {'4: spot linker'};...
        {'8: labelgui - 1'};...
        {'16: tag tracker'};...
        {'32: labelgui - 2'}];
amgHandles.job(activeJobNum).dataProperties.help.cellCycle = [...
        {'cell cycle'};...
        {'1: G1'};...
        {'2: Start'};...
        {'4: S-Phase'};...
        {'8: G2'};...
        {'16: Metaphase'};...
        {'32: Anaphase'};...
        {'64: Telophase'}];
amgHandles.job(activeJobNum).dataProperties.help.strain = [...
        {'strains'};...
        {'1: WT'};...
        {'2: dam1-1'};...
        {'4: dam1-11'};...
        {'8: ipl1-321'};...
        {'16: ipl1-1'};...
        {'32: ndc10-1'};...
        {'64: ndc80-1'};...
        {'128: bik1'};...
        {'256: bim1'}];
amgHandles.job(activeJobNum).dataProperties.help.drugs = [...
        {'drugs'};...
        {'1: nocodazole'};...
        {'2: benomyl'};...
        {'4: hydroxyuracyl'};...
        {'8: alpha-factor'}];

%write editPropertiesGUIDefaults
amgHandles.epguiDefaults.movieProperties.strain =  strainCell;
amgHandles.epguiDefaults.movieProperties.cellCycle = ccCell;
amgHandles.epguiDefaults.movieProperties.drugs = drugCell;
amgHandles.epguiDefaults.movieProperties.tempValue = get(handles.edit_temperature_PD,'Value');
amgHandles.epguiDefaults.background.hList1 = handles.correctBackground1;
amgHandles.epguiDefaults.background.hList2 = handles.correctBackground2;
pdList1 = get(handles.backgroundH(3),'String');
pdList2 = get(handles.backgroundH(4),'String');
%add the "select movie #x" if necessary
if ~strcmp(pdList1{1},'select movie #1')
    pdList1 = [{'select movie #1'};pdList1];
end
if ~strcmp(pdList2{1},'select movie #2')
    pdList2 = [{'select movie #1'};pdList2];
end
amgHandles.epguiDefaults.background.pdList1 = pdList1;
amgHandles.epguiDefaults.background.pdList2 = pdList2;

amgHandles.epguiDefaults.background.framesS = get(handles.backgroundH(6),'String');
amgHandles.epguiDefaults.background.framesL = get(handles.backgroundH(7),'String');

%write dataProperties to file
dataProperties = amgHandles.job(activeJobNum).dataProperties;
save('tmpDataProperties','dataProperties');


%save job data
guidata(amgHandles.AMG,amgHandles);
