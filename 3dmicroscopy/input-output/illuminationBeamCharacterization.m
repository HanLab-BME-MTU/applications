%A script for iterative alignment of the laser illumination and the
%variable slit for adjustment of the excitation numerical aperture.  The
%camera is intended to be a ImagingSource DMK21BU04 CCD.  If the beam is
%rotated clockwise in the propagation direction, then you need to move the
%slit aperture stage counter clockwise.  Very fine movements are sufficient
%to have a significant beam impact on alignment.
% AO 0 = LaserPower.  
% AO 1 = Remote Piezo.
% Y Piezo Scan Settings - 59.838 microns/volt.
% 0 offset is 5 volts
% Kevin Dean, 2015-12-03.

%% Camera registration and initiation.  Camera is 8-bit, and has a 5.6 micron pixel size.  USB 2.0.
clear all
clc
imaqregister('C:\Program Files (x86)\TIS IMAQ for MATLAB R2013b\x64\TISImaq_R2013.dll');
vid = videoinput('tisimaq_r2013', 1, 'Y800 (640x480)');
src = getselectedsource(vid);
vid.FramesPerTrigger = 1;

%% Initiate National Instruments Session.  Depending on computer, may be Dev1 or Dev2.
% Positions y-piezo at 5V and triggers the laser TTL high.

devices=daq.getDevices();  

s = daq.createSession('ni'); 

s.Rate=1E6;

addAnalogOutputChannel(s,'Dev2',0:1,'Voltage'); 

queueOutputData(s,[5 5]);

s.startForeground();


%% Optimize Laser Power.  Only use with AOTF (never go above 5V).  Otherwise input LaserPower=5;

laserPowerOutput=0:0.05:.1; 

outputSize=length(laserPowerOutput);

remotePiezoOutput=ones(outputSize,1)*5; 

clear Stack

imageIntensity=zeros(outputSize,1);

for i=1:outputSize
    queueOutputData(s,[laserPowerOutput(:,i) remotePiezoOutput(i,:)]);
    s.startForeground();
    
    if i == 1
        pause(0.1)
    end
    
    temp=double(getsnapshot(vid));
    
    imageIntensity(i)=max(max(temp));
    
    pause(0.1);
end

i = imageIntensity > 250;

imageIntensity(i)=[];

[int loc] = max(imageIntensity);

LaserPower=laserPowerOutput(loc)

clear i imageIntensity int laserPowerOutput loc outputSize remotePiezoOutput

%% Measure Beam Profile
zStep = 0.01;

remotePiezoOutput=[4.4:zStep:5.7]; %4.8 to 5.2

outputSize=length(remotePiezoOutput);

laserPowerOutput=ones(outputSize,1)*LaserPower;

clear Stack

Stack=zeros(max(size(remotePiezoOutput)),480,640);

for i=1:max(size(remotePiezoOutput))
    queueOutputData(s,[laserPowerOutput(i,:) remotePiezoOutput(:,i)]);
    
    s.startForeground();
        
    if i == 1
        pause(0.1)
    end
    
    temp=double(getsnapshot(vid));
    
    Stack(i,:,:)=temp;
    
    pause(0.1);
end

for i = 1:size(Stack,1)
    imageMax(i) = max(max(Stack(i,:,:)));
end

[imIntensity inFocus] = max(imageMax);


% Identify location of beam and plot the axial cross-section.

temp=squeeze(Stack(ceil(outputSize/2),:,:));

beamProfile=sum(temp,1);  %2 if horizontal, 1 if vertical.

[int loc] = max(beamProfile);

%clear temp beamProfile int

temp=squeeze(Stack(:,loc,:));

temp(:)=0;

for i=1:30; %Mathematical dither of beam in lateral dimension.
    temp=temp+squeeze(Stack(:,loc+i,:)); 
end


% Measure beam FOV FWHM

beamProfile=temp(:,loc);

beamProfileMin=min(beamProfile);

beamProfile=beamProfile-beamProfileMin;

[xData, yData] = prepareCurveData([], beamProfile );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];
%opts.StartPoint = [max(yData) 20 10];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts )

% Plot fit with data.
figure( 'Name', 'Axial Beam Cross-Section' );
subplot(2,2,2);
plot( fitresult, xData, yData );
legend('beamProfile', 'Gaussian Fit', 'Location', 'NorthEast' );
% Label axes
ylabel Intensity
grid on

subplot(2,2,1);
imshow(squeeze(Stack(inFocus,:,:)),[])
colormap jet

subplot(2,2,3:4);
imshow(temp,[]);
colormap jet

FieldOfView=fitresult.c1*2.35*zStep*60.838/sqrt(2)

%%

%axialProfile=squeeze(temp(floor(size(temp,1)/2),:));
axialProfile=squeeze(temp(61,:));

[xData, yData] = prepareCurveData( [], axialProfile-min(axialProfile) );

% Set up fittype and options.
ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts )

% Plot fit with data.
figure( 'Name', 'Axial Beam Cross-Section' );
plot( fitresult, xData, yData );
legend('beamProfile', 'Gaussian Fit', 'Location', 'NorthEast' );
% Label axes
ylabel Intensity
grid on

AxialResolution=fitresult.c1*2.35*.16/sqrt(2)

