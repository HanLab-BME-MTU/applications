function fsmVectorAnalysis(imgNum,roi,displ,scale,d0,useDiv,output,displROI,invertROI)
% fsmVectorAnalysis reports several statistics on vector fields (obtained from FSM)
%
% SYNOPSIS      fsmVectorAnalysis(imgNum,roi,displ,scale,d0,useDiv,output,displROI,invertROI)
%
% REMARK        Run fsmVectorAnalysis through the graphical user interface fsmCenter 
%               ------------------------------------------------------------------
%
% INPUT                     imgNum   :number of images to be analyzed 
%                           roi      : [ 0|1 0|1 0|1]  
%                           roi(1)   : if 1, allows the user to draw a region of interest
%                           roi(2)   : if 1, allows the user to load a saved region of interest
%                           roi(3)   : if 1, saves the drawn region of interest
%                                      (only if roi(1) is 1)
%                           If roi(1) is 1, roi(2) should be 0 and vice versa.
%               displ     : [ 0|1 0|1 0|1 0|1 0|1 ] turns on and off some display
%                           displ(1) : if 1, displays raw vector field
%                           displ(2) : if 1, displays interpolated vector field
%                           displ(3) : if 1, displays noise vectors
%                           displ(4) : if 1, displays map of goodness for raw vs.
%                                      interpolated field match and SNR map
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

if nargin~=9
    error('9 input parameters expected');
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
    INTERP_CALC=1;    % To calculate noise vectors, interpolated vectors are needed
    NOISE_CALC=1;
    NOISE_DISPLAY=1;
else
    NOISE_CALC=0;
    NOISE_DISPLAY=0;
end

if displ(4)==1        % Noise analysis - both interpolated and noise vectors are needed
    INTERP_CALC=1;
    NOISE_CALC=1;
    ERROR_CALC=1;
else
    ERROR_CALC=0;
end

if ~any([RAW_DISPLAY INTERP_CALC NOISE_CALC ERROR_CALC])
    uiwait(msgbox('Nothing to display.','help','modal'));
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SELECT IMAGE AND EXPERIMENT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select image
[fName,dirName] = uigetfile(...
    {'*.tif;*.tiff;*.jpg;*.jpeg','Image Files (*.tif,*.tiff,*.jpg,*.jpeg)';
    '*.tif','TIF files (*.tif)'
    '*.tiff','TIFF files (*.tiff)'
    '*.jpg;','JPG files (*.jpg)'
    '*.jpeg;','JPEG files (*.jpeg)'
    '*.*','All Files (*.*)'},...
    'Select image');
if(isa(fName,'char') & isa(dirName,'char'))
    img=nrm(imread([dirName,fName]),1);
    % Store image size
    imgSize=size(img);
else
    return 
end
% Read frame number
[path,body,no,ext]=getFilenameBody([dirName,fName]);
selectedFrame=str2num(no);

% Select MPM
[mpmFileName,mpmDirName] = uigetfile(...
    {'mpm.mat','mpm.mat'},...
    'Select mpm.mat');
if(isa(mpmFileName,'char') & isa(mpmDirName,'char'))
   load([mpmDirName,mpmFileName]);
else
   return 
end

%
% Calculate which is the frame to be read out of M relative to the first image analyzed
%

% Create structure for filenames
file.name='';

% Open file
fid=fopen([mpmDirName,filesep,'parameters.txt']);
c=0;
while not(feof(fid))
 
    % Read new line
    tline=fgetl(fid);
    
    % Break if end of file
    if ~ischar(tline)
        break
    end
    
    if isempty(tline)
        % Jump to next line
        continue
    end

    if findstr(tline,'First image name and path :')
        c=c+1;
        file(c).name=tline(29:end);
    end
    
end

% Close file
fclose(fid);

% Calculate frame
if c==0
    error('No valid file name for the first image has been found in the file ''parameters.txt''.');
end
fileName=file(c).name;
[path,body,no,ext]=getFilenameBody(fileName);

% Load file names and check input parameter imgNum
if imgNum>1
    outFileList=getFileStackNames([dirName,fName]);
    len=length(outFileList);
    if imgNum>len
        imgNum=len;
    end
end

% Select correct frame in MPM to read
frame=selectedFrame-str2num(no)+1;
if frame<0 | frame>size(M,3)
   error('The selected image has not been analyzed.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Start analysis
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set constants
loadROIAgain=1;   % If a roi is loaded, one cand decied whether to use the same for ALL images, or 
                  % to load one for each new frame - the user will be asked
askLoadAgain=1;   % Ask user?
drawROIAgain=1;   % Draw the roi for each frame, or only once and the use it for all images - the
                  % user will be asked
askDrawAgain=1;   % Ask user?

for c1=1:imgNum
    
    currFrame=frame+c1-1;
    if imgNum>1
        % Load current image
        img=imread(char(outFileList(c1)));
    end
    
    % Extract corresponding frame
    Mm=M(:,:,currFrame);
    
    % Extract vectors
    Mv=Mm(find(Mm(:,1)~=0 & Mm(:,3)~=0),:);
    
    % Sort Mv
    Mv=sortrows(Mv,1:2);
    
    % Set figure handle to zero
    h=0;
    
    if roi(1)==1 & roi(2)==0
        
        if drawROIAgain==1
            
            % Show image
            h=figure; imshow(img,[]);
            
            % Add some explanations...
            set(h,'NumberTitle','off','Name','Please draw your Region Of Interest...');
            
            % Draw roi
            try
                [bw,x,y]=roipoly;
            catch
                uiwait(msgbox('No polygon was selected. Quitting.','Error','modal'));
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
    
            if askDrawAgain==1 & imgNum>1
                button = questdlg('Do you want to use this ROI for all images?',...
                    'User input requested','Yes','No','Yes');
                if strcmp(button,'Yes')
                    drawROIAgain=0;
                elseif strcmp(button,'No')
                    drawROIAgain=1;
                end
                askDrawAgain=0;
            end
            
        end
        
    elseif roi(1)==0 & roi(2)==1
        
        if loadROIAgain==1
            
            % Load roi
            [roiFileName,roiDirName] = uigetfile(...
                {'*.mat','*.mat'},...
                'Load roi...');
            if(isa(roiFileName,'char') & isa(roiDirName,'char'))
                load([roiDirName,roiFileName]);
                if ~(exist('y') & exist('x'))
                    error('The loaded roi is not valid');
                end
            else
                return 
            end
        
            if askLoadAgain==1  & imgNum>1
                button = questdlg('Do you want to use this ROI for all images?',...
                    'User input requested','Yes','No','Yes');
                if strcmp(button,'Yes')
                    loadROIAgain=0;
                elseif strcmp(button,'No')
                    loadROIAgain=1;
                end
                askLoadAgain=0;
            end
        end
        
    elseif roi(1)==1 & roi(2)==1
        error('The region of interest can be either drawn or loaded, not both');
        
    else
        % Roi not drawn and not loaded
        y=[];x=[];
    end
    
    % Cut vectors not belonging to the roi
    if ~isempty(y)
        try
            index=inpolygon(Mv(:,1),Mv(:,2),y,x);
        catch
            % Restart fsmCenter if it exists
            auxH=findobj(0,'Type','figure','Tag','fsmCenter');
            close(auxH);
            fsmCenter;
            uiwait(msgbox('No polygon selected.','Warning','warn'));
            return
        end
        if invertROI==1
            Mv(index,:)=[];        
        else
            Mv(~index,:)=[];
        end
        % Set flag
        POLYGON_DISPLAY=(1 & displROI); % If the user decided not to display the polygon
    else
        POLYGON_DISPLAY=0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % ANALYSIS AND PLOTS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Display image if needed
    if displ(5)==1 & (any([RAW_DISPLAY INTERP_DISPLAY NOISE_DISPLAY])==1)
        
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
    % CALCULATE AND PLOT INTERPOLATED VECTORS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if INTERP_CALC==1
        
        % Extract interpolation points
        Iyx=unique(Mv(:,1:2),'rows');
        
        % Interpolate - get the deterministic part of the signal
        if useDiv==1
            [divM,d0]=vectorFieldDiv(Mv,Iyx,d0_init,[]);
            d0=updateD0FromDiv(divM,d0,1,size(Iyx,1),size(Iyx,1));
        end
        Md=vectorFieldInterp(Mv,Iyx,d0,[]);
        
        % Sort Md
        Md=sortrows(Md,1:2);
        
        if INTERP_DISPLAY==1
            
            % Plot interpolated field
            h=vectorFieldPlot(Md,h,[],scale);
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % CALCULATE AND PLOT NOISE VECTORS
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if NOISE_CALC==1
        
        % "Calculate" noise vector as the vectorial difference between the original vectors and the
        % interpolated ones (deterministic part)
        Ms=zeros(size(Mv));
        Ms(:,1:2)=Md(:,3:4);
        Ms(:,3:4)=Mv(:,3:4);
        
        if scale~=1
            
            Mss=Ms;
            Mss(:,1:2)=Md(:,1:2)+scale*(Md(:,3:4)-Md(:,1:2)); % The noise vector starts at the end of the the interpolated
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
                str=[str,', ',char(colors(2)),' -> interpolated']; % Color 2 used for interpolated
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(3)),' -> noise']; % Color 3 used for interpolated
                end
            else
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(2)),' -> noise']; % Color 2 used for interpolated
                end
            end
        else
            if INTERP_DISPLAY==1
                str=[str,', ',char(colors(1)),' -> interpolated'];
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(2)),' -> noise']; % Color 2 used for interpolated
                end
            else
                if NOISE_DISPLAY==1
                    str=[str,', ',char(colors(1)),' -> noise']; % Color 1 used for interpolated
                end
                
            end      
        end
        title(str);
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
        % Match raw vs. interpolated vectors
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
            title('Match raw vs. interpolated vectors - Red: good match, blue: bad match');
            
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
            title('Match raw vs. interpolated vectors - smaller circles, better match');
            
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
            
        end
    end
    
    if ERROR_CALC==1
        if imgNum==1
            fprintf(1,'File name                        : %s\n',[dirName,fName]);
        else
            fprintf(1,'File name                        : %s\n',char(outFileList(c1)));
        end
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
