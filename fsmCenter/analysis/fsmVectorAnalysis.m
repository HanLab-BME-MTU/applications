function fsmVectorAnalysis(projDir,nAvg,roi,displ,scale,d0,useDiv,output,displROI,invertROI)
% fsmVectorAnalysis reports several statistics on vector fields (obtained from FSM)
%
% SYNOPSIS      fsmVectorAnalysis(projDir,nImages,roi,displ,scale,d0,useDiv,output,displROI,invertROI)
%
% REMARK        Run fsmVectorAnalysis through the graphical user interface fsmCenter 
%               ------------------------------------------------------------------
%
% INPUT         projDir   : string pointing to the project directory (where the fsmParam.mat file is 
%                           located). Pass projDir=[] to manually select fsmParam.mat via a dialog.
%               nImages   : number of images to be analyzed 
%               nAvg      : number of frames for time averaging (must be ODD)
%               roi       :  [ 0|1 0|1 0|1]  
%                           roi(1)   : if 1, allows the user to draw a region of interest
%                           roi(2)   : if 1, allows the user to load a saved region of interest
%                           roi(3)   : if 1, saves the drawn region of interest
%                                      (only if roi(1) is 1)
%                           If roi(1) is 1, roi(2) should be 0 and vice versa.
%               displ     : [ 0|1 0|1 0|1 0|1 0|1 ] turns on and off some display
%                           displ(1) : if 1, displays raw vector field
%                           displ(2) : if 1, displays averaged vector field
%                           displ(3) : if 1, displays noise vectors
%                           displ(4) : if 1, displays map of goodness for raw vs.
%                                      averaged field match and SNR map
%                           displ(5) : if 1, displays image (with all vector fields overlaid)
%                scale    : scaling factor for vector field display
%                d0       : correlation length of the interpolator
%                useDiv   : [ 0|1 ]  : if 1, uses the divergence field to locally adapt
%                           the correlation length d0 for the interpolator
%                output   : [ 1|2 ]  : if 1, displays color-coded error maps, if 2 with
%                           circles with radius proportional to the magnitude of the errors
%                displROI : [ 0|1 ]  : if 1, overlays the drawn polygon
%                invertROI: [ 0|1 ]  : if 1, it selects all vectors OUTSIDE the ROI,
%                                      instead of INSIDE.
%                
% OUTPUT         None  
%
% DEPENDENCES   fsmWaveAnalysis uses { }
%               fsmMain is used by { fsmCenter }
%
% Aaron Ponti, May 15th, 2003

global uFirst uLast

if nargin~=10
    error('10 input parameters expected');
end

if mod(nAvg,2)==0
    error('The number of frames ''nAvg'' must be ODD.')
end

% Store current directory
oldDir=cd;

% Store initial d0
d0_init=d0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SET FLAGS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if displ(1)==1
   RAW_DISPLAY=1; 
else
   RAW_DISPLAY=0; 
end

if displ(2)==1
    INTERP_CALC=1;
    INTERP_DISPLAY=1;
else
    INTERP_CALC=0;
    INTERP_DISPLAY=0;
end

if displ(3)==1
    INTERP_CALC=1;    % To calculate noise vectors, averaged vectors are needed
    NOISE_CALC=1;
    NOISE_DISPLAY=1;
else
    NOISE_CALC=0;
    NOISE_DISPLAY=0;
end

if displ(4)==1        % Noise analysis - both averaged and noise vectors are needed
    INTERP_CALC=1;
    NOISE_CALC=1;
    ERROR_CALC=1;
else
    ERROR_CALC=0;
end

if ~any([RAW_DISPLAY INTERP_CALC NOISE_CALC ERROR_CALC])
    errordlg('Nothing to display...','Error','modal');
    return
end

POLYGON_DISPLAY=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOOK FOR AND CHECK fsmParam.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fsmParam,status]=fsmPostGetFsmParam(projDir);
if status==0
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK THAT THE TRACKING MODULE IS UP-TO-DATE AND LOAD mpm.mat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(fsmParam.track,'uptodate')

    % Check that the tracking module is up-to-date
    if fsmParam.track.uptodate==0
        
        % The tracking module is not up-to-date. Inform the user and return    
        errordlg('The file ''mpm.mat'' is not up-to-date. Please re-run the Tracking module in SpeckTackle.','Error','modal');
        return

    else
        uptodate=1;    
    end
    
else
    
    % Old version of fsmParam.mat. Inform the user that he/she has to make sure that everything is up-to-date
    uiwait(msgbox('Since ''fsmParam.mat'' has been created by an old version of SpeckTackle, I cannot make sure that ''mpm.mat'' is up-to-date. Continue at your own risk.','Warning'));
    uptodate=-1;
    
end

% Load mpm.mat
if exist([projDir,filesep,'mpm.mat'],'file');
    load([projDir,filesep,'mpm.mat']);
    clear MPM;
else
    errordlg('Could not find ''mpm.mat''.','Error','modal');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GET THE LIST OF IMAGE FILE NAMES AND FIRST/LAST INDICES FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imageFileList,imageFirstIndex,imageLastIndex,status]=fsmPostGetImageListFromFsmParam(fsmParam);
if status==0
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHECK THAT THE SIZE OF M MATCHES THE NUMBER OF ANALYZED IMAGES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nImages=imageLastIndex-imageFirstIndex+1;
if size(M,3)+1~=nImages
    if uptodate==-1
        errordlg('The tracking module is not up-to-date. Please re-run the tracking module in SpeckTackle.','Error','modal');
    else
        errordlg('Even though the tracking module appears to be up-to-date, it is NOT. THIS IS A BUG. Please REPORT it.','Error','modal');
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ASK THE USER TO SPECIFY THE RANGE OF IMAGES HE/SHE WANTS TO ANALYZE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The minimum number of frames to be selected depends on the value assigned to nAvg by the user
if nAvg==0
    minN=1; % This causes the all stack to be averaged into one speed map
else
    minN=nAvg-1;
end

guiH=fsmTrackSelectFramesGUI; ch=get(guiH,'Children');
set(findobj('Tag','pushOkay'),'UserData',minN); % At least n-1 frames must be considered
titleDlg='Select frame pairs to be processed:';
set(findobj('Tag','editFirstFrame'),'String',num2str(1));
set(findobj('Tag','editLastFrame'),'String',num2str(nImages-1));
set(findobj('Tag','SelectFramesGUI'),'Name',titleDlg);
sSteps=[1/((nImages-1)-1) 1/((nImages-1)-1)];
set(findobj('Tag','sliderFirstFrame'),'SliderStep',sSteps,'Max',nImages-1,'Min',1,'Value',1);
set(findobj('Tag','sliderLastFrame'),'SliderStep',sSteps,'Max',nImages-1,'Min',1,'Value',nImages-1);
waitfor(guiH); % The function waits for the dialog to close (and to return values for uFirst and uLast)

if uFirst==-1
    return % The user closed the dialog
end

% Update n depending on the selection
if nAvg==0
    nAvg=uLast-uFirst+1;
end

if nAvg>uLast
    nAvg=uLast;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CREATE NUMER AND LIST OF IMAGES, AND FRAMES TO BE PROCESSED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Index of the first and last matrices in M to be analyzed
firstMatrix=uFirst;
lastMatrix=uLast;

% Calculate the corresponding indices for imageFileList and for the outputs
imageIndices=[firstMatrix:lastMatrix-nAvg+1]+fix(nAvg/2);

% Update number of images to be processed
nImages=length(imageIndices);

% Update imageFileList
imageFileList=imageFileList(imageIndices,:);

% Crop the frames to be processed from M
M=M(:,:,firstMatrix:lastMatrix);

% Load first image
img=imread(char(imageFileList(1,:)));

% Store image size
imgSize=size(img);

% If needed, ask the user to specify an output directory
if nImages>1
    % Select output dir
    outputdir=uigetdir('','Select directory to save vector maps to.');
    if outputdir==0 % The user clicked on cancel
        disp('Aborted by the user.');
        return
    end
    
    % String format
    [path,outputFileName,no,ext]=getFilenameBody(imageFileList(1,:));
    s=length(no);
    strg=sprintf('%%.%dd',s);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DRAW | LOAD, SAVE A POLYGON IF NEEDED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if roi(1)==1 & roi(2)==0
    
    % Show image
    h=figure; imshow(img,[]);
    
    % Add some explanations...
    set(h,'NumberTitle','off','Name','Please draw your Region Of Interest...');
    
    % Draw roi
    try
        [bw,x,y]=roipoly;
    catch
        errordlg('No polygon was selected. Quitting.','Error','modal');
        return
    end
    
    % Set coordinates outside of the image to the image borders + (or -) 1
    x(find(x<1))=0; x(find(x>imgSize(2)))=imgSize(2)+1;
    y(find(y<1))=0; y(find(y>imgSize(1)))=imgSize(1)+1;
    
    % Close figure
    close(h);
    
    % Save roi
    if roi(3)==1
        
        % Specify path and filename
        [saveFileName,savePath]=uiputfile(...
            {'*.mat','*.mat'},...
            'Save polygon as...');
        if(isa(saveFileName,'char') & isa(savePath,'char'))
            if isempty(findstr('.',saveFileName))
                saveFileName=[saveFileName,'.mat'];
            else
                if ~strcmp(saveFileName(end-3:end),'.mat')
                    saveFileName(end-3:end)='.mat';
                end
            end
            eval(['save ',[savePath,filesep,saveFileName],' y x']);
        else
            disp('ROI not saved.');
        end
        
    end
    
elseif roi(1)==0 & roi(2)==1
    
    % Load roi
    [roiFileName,roiDirName] = uigetfile(...
        {'*.mat','*.mat'},...
        'Load roi...');
    if(isa(roiFileName,'char') & isa(roiDirName,'char'))
        load([roiDirName,roiFileName]);
        if ~(exist('y') & exist('x'))
            errordlg('The loaded roi is not valid.','Error','modal');
            return
        end
    else
        return 
    end
    
elseif roi(1)==1 & roi(2)==1
    error('The region of interest can be either drawn or loaded, not both');
    
else
    % Roi not drawn and not loaded
    y=[];x=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% APPLY POLYGON IF REQUESTED BY THE USER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(y) & ~isempty(x)
    
    if displROI==1
        PPOLYGON_DISPLAY=1;
    end
    
    for c1=1:size(M,3)
        
        % Extract corresponding frame
        Mm=M(:,:,c1);
        
        % Extract vectors
        Mv=Mm(find(Mm(:,1)~=0 & Mm(:,3)~=0),:);
        
        % Cut vectors not belonging to the roi
        index=inpolygon(Mv(:,1),Mv(:,2),y,x);
        
        % If the user selected 'Invert ROI' kill the vectors WITHIN the polygon
        if invertROI==1
            Mv(index,:)=[];        
        else
            Mv(~index,:)=[];
        end
        
        % Store the selected vectors
        M(:,:,c1)=0;
        M(1:size(Mv,1),1:4,c1)=Mv;
        
    end
    
    % Resize M
    cM=M>0; [i j k]=find(cM); M=M(1:max(i),:,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INTERPOLATE VECTOR FIELDS IF NEEDED
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if INTERP_CALC==1
    
    % Initialize empty Mt stack where to store the averaged vector fields
    emptyMt=[0 0 0 0];
    Mt=emptyMt;
    
    currentFrame=0;
    
    % Interpolation frame shift
    intFrame=fix(nAvg/2)+1;
        
    for c1=1:nImages
        
        currentM=[];
        
        % To average over time, append raw vectors from the nAvg frames
        for c2=c1:c1+nAvg-1
            
            % Extract corresponding frame
            Mm=M(:,:,c2);
        
            % Extract vectors
            Mv=Mm(find(Mm(:,1)~=0 & Mm(:,3)~=0),:);

            %Append
            currentM=cat(1,currentM,Mv);    
            
        end
        
        % The interpolation points are the coordinates of the central frame
        Mm=M(:,:,c1+intFrame-1);
        Iyx=Mm(find(Mm(:,1)~=0 & Mm(:,3)~=0),1:2);
            
        % Interpolate - get the deterministic part of the signal
        try
            if useDiv==1
                [divM,d0]=vectorFieldDiv(currentM,Iyx,d0_init,[]);
                d0=updateD0FromDiv(divM,d0,1,size(Iyx,1),size(Iyx,1));
            end
            Md=vectorFieldInterp(currentM,Iyx,d0,[]);
        catch
            errordlg('Sorry, OUT OF MEMORY. Either reduce the number of frames for time averaging or analyze only a subset of the vector field (by drawing a ROI).','Error','modal');
            return
        end
        
        % Store Md
        currentFrame=currentFrame+1;
        Mt(1:size(emptyMt,1),1:size(emptyMt,2),currentFrame)=emptyMt;
        Mt(1:size(Md,1),1:size(Md,2),currentFrame)=Md;
            
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ANALYSIS AND PLOTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c1=1:nImages %firstMatrix:lastMatrix
    
    % Current index for raw data (taking into account the number of frames for averaging)
    currentRawFrame=c1+fix(nAvg/2);
    
    % Extract corresponding frame
    Mm=M(:,:,currentRawFrame);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CURRENT RAW VECTORS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mv=Mm(find(Mm(:,1)~=0 & Mm(:,3)~=0),:); 
    
    % Set figure handle to zero
    h=0;
    
    % Display image if needed
    if displ(5)==1 & (any([RAW_DISPLAY INTERP_DISPLAY NOISE_DISPLAY])==1)
        
        % Load current image
        img=imread(char(imageFileList(c1,:)));
        
        % Show image
        h=figure; imshow(img,[]);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % PLOT RAW VECTORS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if RAW_DISPLAY==1
        
        % Plot vector field
        h=vectorFieldPlot(Mv,h,[],scale); % h is zero if no figure exists, othersiwe it is equal to its handle
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % PLOT INTERPOLATED VECTORS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if INTERP_CALC==1 & INTERP_DISPLAY==1
        
        % Extract corresponding frame
        Mm=Mt(:,:,c1);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CURRENT INTERPOLATED VECTORS
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Md=Mm(find(Mm(:,1)~=0 & Mm(:,3)~=0),:);
        
        % Plot averaged field
        h=vectorFieldPlot(Md,h,[],scale);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CALCULATE AND PLOT NOISE VECTORS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if NOISE_CALC==1
        
        % "Calculate" noise vector as the vectorial difference between the original vectors and the
        % averaged ones (deterministic part)
        Ms=zeros(size(Mv));
        Ms(:,1:2)=Md(:,3:4);
        Ms(:,3:4)=Mv(:,3:4);
        
        if scale~=1
            
            Mss=Ms;
            Mss(:,1:2)=Md(:,1:2)+scale*(Md(:,3:4)-Md(:,1:2)); % The noise vector starts at the end of the the averaged
            Mss(:,3:4)=Mv(:,1:2)+scale*(Mv(:,3:4)-Mv(:,1:2)); % (deterministic) vector and ends at the end of the raw vector 
            
        else
            
            Mss=Ms;
            
        end
        
        if NOISE_DISPLAY==1
            
            % Plot stochastic field (noise)
            h=vectorFieldPlot(Mss,h,[],1);
            
        end
        
    end
    
    if any([RAW_DISPLAY INTERP_DISPLAY NOISE_DISPLAY])==1
        
        % Set axes
        axis([1 imgSize(2) 1 imgSize(1)]);
        
        % Draw polygon
        if POLYGON_DISPLAY==1
            hold on;
            plot(x,y,'r-');
        end
        
        % Add title
        if displ(5)==1
            colors={'yellow','red','blue'};
        else
            colors={'black','red','blue'};
        end
        
        str=['Scaled ',num2str(scale),'x'];
        if RAW_DISPLAY==1
            str=[str,', ',char(colors(1)),' -> raw']; % Color 1 used for raw
            if INTERP_DISPLAY==1
                str=[str,', ',char(colors(2)),' -> averaged']; % Color 2 used for averaged
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(3)),' -> noise']; % Color 3 used for averaged
                end
            else
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(2)),' -> noise']; % Color 2 used for averaged
                end
            end
        else
            if INTERP_DISPLAY==1
                str=[str,', ',char(colors(1)),' -> averaged'];
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(2)),' -> noise']; % Color 2 used for averaged
                end
            else
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(1)),' -> noise']; % Color 1 used for averaged
                end
                
            end      
        end
        title(str);
        
        % Save the image to disk if needed
        if nImages>1
            indxStr=sprintf(strg,imageIndices(c1));
            fname=[outputdir,filesep,'flowMap_d0=',num2str(d0_init),'_frames=',num2str(nAvg),'_',indxStr,'.tif'];
            print(gcf,'-dtiffnocompression',fname);
            
            % Close image
            close(h);
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % ERROR CALCULATIONS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ERROR_CALC==1
        
        % Extract vectors
        d=Md(:,3:4)-Md(:,1:2); % Interpolated vector field
        v=Mv(:,3:4)-Mv(:,1:2); % Raw vector field
        s=Ms(:,3:4)-Ms(:,1:2); % Noise field
        
        % Vector lengths
        ld=sqrt(d(:,1).^2+d(:,2).^2);
        lv=sqrt(v(:,1).^2+v(:,2).^2);
        ls=sqrt(s(:,1).^2+s(:,2).^2);
        
        %
        % Match raw vs. averaged vectors
        %
        
        % Error on lengths
        eL=abs(lv-ld)./ld;                      % Relative error on length
        
        % Angles - vectorize
        dy=d(:,1); dx=d(:,2); vy=v(:,1); vx=v(:,2);
        num=sqrt(vy.^2+vx.^2).*sqrt(dy.^2+dx.^2);
        num(find(num==0))=eps;
        eA=acos((vy.*dy+vx.*dx)./num)/pi;
        
        % Total error
        eT=sqrt(eL.^2+eA.^2);
        
        %
        % Calculate signal to noise ratio
        %
        
        warning off
        SNR=ld./ls;
        warning on
        
        % Construct matrices for display
        E=[Mv(:,1:2) eT];
        S=[Mv(:,1:2) SNR];
        
        %
        % Display errors
        %
        
        if output==1
            
            %
            % Colormap
            %
            
            % Open new figure for vector matches
            h=figure; 
            
            % Create empty map
            match=zeros(imgSize);
            
            % Put errors
            for i=1:size(E,1)
                match(E(i,1),E(i,2))=31/E(i,3);
            end
            
            % Calculate the approximate mean error distance and sigma for the gaussian kernel
            sigma=sqrt(prod(imgSize)/size(E,1))/3; % Doesn't take into accout whether the error 
            % distribution occupies the whole image surface or not
            
            % Filter result with a gaussian kernel and the calculated sigma
            match=Gauss2D(match,sigma);
            
            % Apply color map
            img3C_map=applyColorMap(img,match,[0 30],jet(31),16);
            
            % Show image
            imshow(img3C_map,[]);
            title('Match raw vs. averaged vectors - Red: good match, blue: bad match');
            
            % Save the image to disk if needed
            if nImages>1
                indxStr=sprintf(strg,imageIndices(c1));
                fname=[outputdir,filesep,'matchMap_d0=',num2str(d0_init),'_frames=',num2str(nAvg),'_',indxStr,'.tif'];
                print(gcf,'-dtiffnocompression',fname);
                
                % Close image
                close(h);

            end

            % Open new figure for SNR
            h=figure; 
            
            % Create empty map
            SNRmap=zeros(imgSize);
            
            % Put errors
            for i=1:size(S,1)
                SNRmap(S(i,1),S(i,2))=31*S(i,3);
            end
            
            % Filter with the same sigma as before
            SNRmap=Gauss2D(SNRmap,sigma);
            
            % Apply color map
            img3C_SNR=applyColorMap(img,SNRmap,[0 30],jet(31),16);
            
            % Show image
            imshow(img3C_SNR,[]);
            title('SNR - Red: good match, blue: bad match');
            
            % Save the image to disk if needed
            if nImages>1
                indxStr=sprintf(strg,imageIndices(c1));
                fname=[outputdir,filesep,'SNRMap_d0=',num2str(d0_init),'_frames=',num2str(nAvg),'_',indxStr,'.tif'];
                print(gcf,'-dtiffnocompression',fname);
                
                % Close image
                close(h);

            end
            
        else
            
            %
            % Circles
            %
            
            % Calculate sin and cos vectors to display circles
            ang=[0:pi/4:2*pi];
            sn_small=sin(ang);
            cs_small=cos(ang);
            ang=[0:pi/10:2*pi];
            sn_med=sin(ang);
            cs_med=cos(ang);
            ang=[0:pi/100:2*pi];
            sn_large=sin(ang);
            cs_large=cos(ang);
            
            % Open new figure for vector matches
            h=figure; 
            
            if displ(5)==1
                
                % Show image
                imshow(img,[]);
                hold on
                
            end
            
            % Plot first error
            if E(1,3)<=2
                plot(E(1,2)+E(1,3)*sn_small,E(1,1)+E(1,3)*cs_small,'r.');
            elseif E(1,3)>2 & E(1,3)<10
                plot(E(1,2)+E(1,3)*sn_med,E(1,1)+E(1,3)*cs_med,'r.');
            else
                plot(E(1,2)+E(1,3)*sn_large,E(1,1)+E(1,3)*cs_large,'r.');
            end            
            hold on;
            
            % Plot the rest
            for i=2:size(E,1)
                if E(i,3)<=2
                    plot(E(i,2)+E(i,3)*sn_small,E(i,1)+E(i,3)*cs_small,'r.');
                elseif E(i,3)>2 & E(i,3)<10
                    plot(E(i,2)+E(i,3)*sn_med,E(i,1)+E(i,3)*cs_med,'r.');
                else
                    plot(E(i,2)+E(i,3)*sn_large,E(i,1)+E(i,3)*cs_large,'r.');
                end            
            end
            axis ij
            title('Match raw vs. averaged vectors - smaller circles, better match');
            
            % Save the image to disk if needed
            if nImages>1
                indxStr=sprintf(strg,imageIndices(c1));
                fname=[outputdir,filesep,'matchMap_d0=',num2str(d0_init),'_frames=',num2str(nAvg),'_',indxStr,'.tif'];
                print(gcf,'-dtiffnocompression',fname);
                
                % Close image
                close(h);

            end
            
            % Open new figure for SNR
            h=figure; 
            
            if displ(5)==1
                
                % Show image
                imshow(img,[]);
                hold on
                
            end
            
            % Plot first error
            if S(1,3)<=2
                plot(S(1,2)+S(1,3)*sn_small,S(1,1)+S(1,3)*cs_small,'r.');
            elseif S(1,3)>2 & S(1,3)<10
                plot(S(1,2)+S(1,3)*sn_med,S(1,1)+S(1,3)*cs_med,'r.');
            else
                plot(S(1,2)+S(1,3)*sn_large,S(1,1)+S(1,3)*cs_large,'r.');
            end            
            hold on;
            
            
            % Plot the rest
            for i=2:size(S,1)
                if S(i,3)<=2
                    plot(S(i,2)+S(i,3)*sn_small,S(i,1)+S(i,3)*cs_small,'r.');
                elseif E(i,3)>2 & E(i,3)<10
                    plot(S(i,2)+S(i,3)*sn_med,S(i,1)+S(i,3)*cs_med,'r.');
                else
                    plot(S(i,2)+S(i,3)*sn_large,S(i,1)+S(i,3)*cs_large,'r.');
                end            
            end
            axis ij
            title('SNR - larger circles, higher SNR');
            
            % Save the image to disk if needed
            if nImages>1
                indxStr=sprintf(strg,imageIndices(c1));
                fname=[outputdir,filesep,'SNRMap_d0=',num2str(d0_init),'_frames=',num2str(nAvg),'_',indxStr,'.tif'];
                print(gcf,'-dtiffnocompression',fname);
                
                % Close image
                close(h);

            end

        end
    end
    
    if ERROR_CALC==1

        fprintf(1,'File name                        : %s\n',char(imageFileList(c1,:)));
        fprintf(1,'Correlation length d0            : %d\n',d0_init);
        fprintf(1,'Number of RAW vectors            : %d\n',size(Mv,1));
        fprintf(1,'Mean RAW vector length           : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(lv),std(lv),100*std(lv)/mean(lv));
        fprintf(1,'Median RAW vector length         : %2.4f\n',median(lv));
        fprintf(1,'Mean INTERPOLATED vector length  : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(ld),std(ld),100*std(ld)/mean(ld));
        fprintf(1,'Median INTERPOLATED vector length: %2.4f\n',median(ld));
        fprintf(1,'Mean NOISE vector length         : %2.4f +/- %2.4f (+/- %2.2f%%)\n',mean(ls),std(ls),100*std(ls)/mean(ls));
        fprintf(1,'Median NOISE vector length       : %2.4f\n',median(ls));
        snr=ld./ls;
        fprintf(1,'Mean / median SNR                : %2.4f +/- %2.4f (+/- %2.2f%%) / %2.4f\n',mean(snr),std(snr),100*std(snr)/mean(snr),median(snr));
        [n,h]=hist(snr,[0.5:max(snr)+0.5]);
        n=cumsum(n); n=n/max(n);
        indx=find(n>=0.95);indx=indx(1);
        fprintf(1,'95%% of all vectors have SNR <= %d\n',indx);
    end
    
end

% Change back to the initial directory
cd(oldDir);
