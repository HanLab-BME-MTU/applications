clear 
clc


% dialog_title = 'Select the directory containing the images to be processed'; % select images
% folder_name = uigetdir('',dialog_title);
folder_name = '/work/bioinformatics/s170480/Data/FA/Cambridge/beads_20/';
addpath(folder_name);

% select current folder
cd(folder_name);

ImagesToRead = dir('*.tif');

%preallocation

ImageCell = cell(1,length(ImagesToRead));

for i=1 : length(ImagesToRead)    
   ImageCell{i} = imread(ImagesToRead(i).name);    
end


optimizer = registration.optimizer.OnePlusOneEvolutionary();% registration.optimizer.RegularStepGradientDescent; % here you can modify the default properties of the optimizer to suit your need/to adjust the parameters of registration.
optimizer.InitialRadius=6.250000e-05;
optimizer.MaximumIterations=100;

[optimizer2, metric]  = imregconfig('monomodal'); % for optical microscopy you need the 'monomodal' configuration.

RegisteredCell = cell(1,length(ImagesToRead));
RegisteredCell{1} = ImageCell{1};
for p = 2:length(ImagesToRead)

   moving = ImageCell{p}; % the image you want to register
   fixed = ImageCell{1}; % the image you are registering with

   % MATLAB
   tform = imregtform(moving,fixed,'rigid',optimizer,metric,'DisplayOptimization',true,'PyramidLevels',5);                   
   tform = affine2d(tform.T);
   RegisteredCell{p} = imwarp(moving, tform,'OutputView', imref2d(size(fixed))); % 
   outputSetMAT(p,:) = tform.T(3,1:2);
   
%  Paper 2008 
   usfac = 20;
   [output, Greg] = dftregistration(fft2(fixed),fft2(moving),usfac);
   RegisteredCell{p} = abs(ifft2(Greg));
   outputSet(p,:) = output;
   
   figure(2);
   imshowpair(fixed, RegisteredCell{p}, 'Scaling','joint');  
     
end



h=figure(33);
pos = get(gcf,'Position');
h   = figure(gcf); hold off;
set(h,'Position',pos);
for k = 1:length(ImagesToRead)
%    img = imread(imgFile{k});
   imshow(RegisteredCell{k},[]);
   set(h,'Position',pos);
   M(k) = getframe;
end


filewritename = fullfile(folder_name, 'newEffp_testNow100.avi');
v = VideoWriter(filewritename);
open(v);
writeVideo(v, M);
close(v);





% outputSet =
% 
%          0         0         0         0
%     0.0930   -0.0000    0.1000   -0.0500
%     0.0932   -0.0000    0.2500    0.0500
%     0.0951   -0.0000    0.4000   -0.1000
%     0.0961    0.0000    0.5000   -0.0500
%     0.0968    0.0000    0.2500   -0.1500
%     0.0971   -0.0000         0    0.2000
%     0.0978    0.0000    0.2000    0.0500
%     0.0973   -0.0000    0.3500   -0.2000
%     0.0992   -0.0000    0.4500    0.2500
%     0.0981   -0.0000    0.4000    0.2500
%     0.0994   -0.0000   -0.0500   -0.0500
%     0.0989    0.0000    0.1500    0.0500
%     0.1007   -0.0000    0.3000    0.3500
%     0.1008   -0.0000   -0.0500    0.1000
% 
% outputSetMAT
% 
% outputSetMAT =
% 
%          0         0
%    -0.2997    0.3424
%     0.3341    0.4059
%    -0.3347    0.5353
%    -0.3513    0.6098
%    -0.3708    0.4311
%     0.3977   -0.2967
%     0.3145    0.3595
%    -0.4138    0.4209
%     0.4185    0.5068
%     0.4319    0.5035
%    -0.3582   -0.2769
%    -0.2835    0.3805
%     0.4969    0.3851
%     0.3658   -0.3409


