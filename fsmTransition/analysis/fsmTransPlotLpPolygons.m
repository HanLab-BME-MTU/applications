function fsmTransPlotLpPolygons(lpPolygons,plotSize)
% fsmTransPlotLpPolygons allows the user to plot polygons in an interactive fashion
%
% fsmTransPlotLpPolygons(lpPolygons,plotSize)
%
% INPUT
%
%          lpPolygons     : polygons generated through the "Create/plot
%                           lamellipodium polygons" button in fsmTransition
%                           Structure:     lpPolygons(1:number of frames).y : y coordinates
%                                                                        .x : x coordinates
%          plotSize       : (optional) size of the original image [y x] 
%
% OUTPUT
%          none
%
% Aaron Ponti, 12/08/2003

% If no plotSize has been passed, calculate it from the coordinates to be plotted
if nargin==1
    mx=0;
    my=0;
    for i=1:length(lpPolygons)
        cx=max([lpPolygons(i).x]);if cx>mx, mx=cx; end;
        cy=max([lpPolygons(i).y]);if cy>my, my=cy; end;
    end
    plotSize=[cy cx];
end

% Set up plot
fH=figure;
set(fH,'Name','Lp polygons');
set(fH,'NumberTitle','off');
axis equal 
axis ij
axis([1 plotSize(2) 1 plotSize(1)]);
cont=1;
n=1:length(lpPolygons);
nm=n(end);
pos=1;
i=n(pos);
acc=0;
marker='r.-';
markerStatus=0;
toPlot=1;

% Go into an interactive loop
while cont==1
    
    %title(['Previous: ''a'' (new plot), ''z'' (hold on),''q'' (no plot); next: ''s'' (new plot), ''x'' (hold on), ''w'' (no plot); ''e'' to exit; ''m'' turn on/off marker; ''c'' clear. Current frame: ',num2str(n(pos))]);
    title(['Current frame: ',num2str(n(pos)),' - Press (h) for help, (e) to exit']);
    try
        t=waitforbuttonpress;
    catch
        cont=0;
        break;
    end
    if t==1
        ch=get(fH,'CurrentCharacter');
        if ch=='h'
            instring={' ', ...
                    '   h:   displays this help', ...
                    ' ', ...
                    'PLOT', ...
                    ' ', ...
                    ' Previous frame:', ...
                    '   a:   resets the axes and plots the previous frame', ...
                    '   z:   plots the previous frame on top of the current', ...        
                    '   q:   moves back to the previous frame without plotting it',...         
                    ' ', ...
                    ' Next frame:', ...
                    '   s:   resets the axes and plots the next frame', ...
                    '   x:   plots the next frame on top of the current', ...        
                    '   w:   moves back to the next frame without plotting it',...
                    ' ', ...
                    'MORE OPTIONS',...
                    ' ', ...
                    '   m:   Turns on and off line markers (and sets the line width of the', ...
                    '           leading edge larger or equal to the transition''s line)',...
                    '   c:   clears the plot'};
            uiwait(msgbox(instring,'Help: keys and their actions'));
            toPlot=0;
        end
        if ch=='m'
            toPlot=1;
            markerStatus=mod(markerStatus+1,2);
            if markerStatus==0
                marker='r.-';
            else
                marker='r';
            end
        end
        if ch=='a' | ch=='z' | ch=='q'
            if ch=='q'
                toPlot=0;
            else
                toPlot=1;
            end
            if ch=='a'
                cla;   % Clear axes
                acc=0; % Turn off accumulation of plots
            elseif ch=='z'
                acc=1; % Turn on accumulation of plots
            else
                % Nothing
            end
            pos=pos-1;
            if pos<1
                pos=nm; % Wrap around
            end
            i=n(pos);
        elseif ch=='s' | ch=='x' | ch=='w'
            if ch=='w'
                toPlot=0;
            else
                toPlot=1;
            end
            if ch=='s'
                cla;   % Clear axes             
                acc=0; % Turn off accumulation of plots
            elseif ch=='x'
                acc=1; % Turn on accumulation of plots
            else
                % Nothing
            end
            pos=pos+1;
            if pos>nm
                pos=1; % Wrap around
            end
            i=n(pos);
        elseif ch=='e'
            close(fH);
            cont=0;
            break % Stop drawing now
        elseif ch=='c'
            cla
            toPlot=0;
        else
        end
        
        % Turns on/off accumulation of plots
        if acc==1
            hold on
        else
            hold off
        end
        
        % If necessary, plots current frame
        if toPlot==1
            plot(lpPolygons(i).x,lpPolygons(i).y,marker);
            axis ij
            axis equal
            axis([1 plotSize(2) 1 plotSize(1)]);
            title(num2str(i));
        end
        
    end
end
