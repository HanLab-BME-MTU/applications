function []=makeCellNetworkMovie(network,target_dir,xLimVal,yLimVal,scale,fps,movieFormat,noErrs,noTicks)
% makeCellNetworkMovie(trackedNet,'cellNetwork',[550 1250],[300 900],1,8,'mov',1,0)
%get the target directory:
if nargin < 2 || isempty(target_dir)
    target_dir = uigetdir('','Please select target directory');
end

if ~isdir(target_dir)
    mkdir(target_dir)
end

if nargin < 6 || isempty(fps)
    fps = 8;
end

movieFileName='mov_cellNetwork';
if nargin >= 7 && (isempty(movieFormat) || strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV'))
    movieFileName = [movieFileName,'.mov'];
    MakeQTMovie('start',movieFileName);
    MakeQTMovie('framerate',fps);
    MakeQTMovie('quality',1);
    movieFormat ='mov';
    doMovie=1;
elseif nargin >= 7 && (strcmp(movieFormat,'avi') == 1 || strcmp(movieFormat,'AVI'))
    movieFileName = [movieFileName,'.avi'];
    movieFormat ='avi';
    doMovie=1;
else
    doMovie=0;
end

if nargin < 8 || isempty(noErrs)
    noErrs = 1;
end

if nargin < 9 || isempty(noTicks)
    noTicks = 0;
end

padZeros=3;
k=1;
for frame=1:length(network)
    % Here starts the image:
    % scrsz = get(0,'ScreenSize');
    % h=figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);
    if ~noErrs || network{frame}.stats.errs==0        
        h=figure();
        if ~isempty(network{frame})
            display(['Plot frame: ',num2str(frame)]);
            
            % plot the cell network:
            plotCellNetwork(network{frame},xLimVal,yLimVal,scale,noErrs,noTicks);
            
            filename = [target_dir,filesep,'cellNetwork_',num2str(frame,['%0.',int2str(padZeros),'d'])];
            % saveas(gcf,[filename,'.tiff'],'tiffn');
            
            print('-depsc2','-loose', 'frame.eps');
            str = ['!convert -density 300 -quiet -colorspace rgb -alpha off -depth 8 frame.eps ',filename,'.tif'];
            eval(str);
            
            % saveas(gcf,[filename, '.eps'], 'psc2');
            % display(['Figure saved to: ',filename,'.tiffn+.eps'])
            
            pause(0.5)
            
            if doMovie==1 && (strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV'))
                MakeQTMovie('addfigure');
            elseif doMovie==1
                M(k) = getframe;
            end
            
            if frame<length(network)
                close(h);
            end
        end
    else
        display('Found at least one node with unresolved force field! To plot this use noErrs=0')
    end
    k=k+1;
end
!rm frame.eps


if doMovie==1 && (strcmp(movieFormat,'mov') == 1 || strcmp(movieFormat,'MOV'))
    MakeQTMovie('finish');
elseif doMovie==1 && (strcmp(movieFormat,'avi') == 1 || strcmp(movieFormat,'AVI'))
    movie2avi(M,movieFileName,'fps',fps);
end