function psf2D=fsm2dPolym(firstTime,psf,polym,nPix,psf2Di,winPos);
% sub function for fsm demo: plots the 2d polymer

polym2D=polym;
psf=psf(ceil(length(psf)/4):ceil(length(psf)*3/4));

% open display window for polymer
if(firstTime)
   polyWinH = figure('Units','points', ...
      'MenuBar','none', ...
      'Name','Polymer', ...
      'NumberTitle','off', ...
      'Resize','on', ...
      'Position',winPos);
	lenPsf=length(psf);
	centPsf=ceil(lenPsf/2);
   
   axHandle = axes('Parent', polyWinH);

	% Create 2D PSF and convolve with polymer
   rd=zeros(lenPsf);
   
   
	%create coordinates
	[rx ry]=ind2sub([lenPsf lenPsf],1:length(rd(:)));
	%compute distance from center
	d = round(sqrt((rx-centPsf).^2+(ry-centPsf).^2));
	for l = 1: length(d)
   	if(d(l)<centPsf)
      	rd(rx(l),ry(l))=psf(centPsf+d(l));
   	end;
	end;
   psf2D = rd / sum(rd(:));
   
cdatamapping = 'scaled';
h= image(1:30,1:100, zeros(100,30), 'BusyAction', 'cancel', ...
   'Parent', axHandle, 'CDataMapping', cdatamapping, ...
   'Interruptible', 'off');
set(axHandle, ...
   'TickDir', 'out', ...
   'XGrid', 'off', ...
   'YGrid', 'off', ...
   'DataAspectRatio', [1 1 1], ...
   'PlotBoxAspectRatioMode', 'auto', ...
   'Visible', 'off');  
%set(axHandle, 'Units', 'normalized', 'Position', [0 0 1 1]);
map = gray(256);
set(polyWinH, 'Colormap', map);
clim = [0 1];
set(axHandle, 'CLim', clim);

else
   polyWinH=findobj(0,'Name','Polymer');
   psf2D=psf2Di;
end;


% 2. Step: convolution
picPolym2D = conv2(polym2D,psf2D);

% generate pixelated picture
pixPicPolym2D=zeros(floor(size(picPolym2D)/nPix));
for(j=1:(size(pixPicPolym2D,2)-1))
   for(i=1:(size(pixPicPolym2D,1)-1))
      pixField=picPolym2D(ceil((i-.5)*nPix+1):floor((i+.5)*nPix),ceil((j-.5)*nPix+1):floor((j+.5)*nPix));
      pixPicPolym2D(i,j) = sum(pixField(:));
   end;
end;

% adjust the mean intensity
maxPicPolym2D = max(pixPicPolym2D(:));
if(maxPicPolym2D > 0)
   pixPicPolym2D = pixPicPolym2D/maxPicPolym2D;
end;

% plot functions
figure(polyWinH);
imH=findobj(polyWinH,'Type','image');
axHandle=findobj(polyWinH,'Type','axes');

xdata = (1:size(pixPicPolym2D,2));
ydata = (1:size(pixPicPolym2D,1));
set(imH,'XData',xdata);
set(imH,'YData',ydata);
set(imH,'CData',pixPicPolym2D);
   