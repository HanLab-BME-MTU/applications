function r3dview(img,spots,mnp)
% R3DVIEW displays *.r3d image series
% does not (yet) support multiple wvls in images

% c: 7/3/01 dT

if(nargin==0)
   [img, fname]=r3dread;
   if isempty(img)
      return;
   end;
else
   fname='img';
end;

% intialize params
iptsetpref('ImshowBorder', 'tight');
img=img/max(img(:));
imSize=size(img);

%check if already an open panel exists
imgFigureH=findobj('Type','figure','Name','View Panel');
contFigH=findobj('Type','figure','Name','r3d viewer');
if ~isempty(contFigH)
    delete(contFigH);
end;
if ~isempty(imgFigureH)
    delete(imgFigureH);
end;

%open control panel
contFigH=r3dviewpanel;

%set first image
imgFigH=uiViewPanelShowImg(img(:,:,1,1,1));
imgH=findobj(imgFigH,'Type','image');
set(imgH,'CData',img(:,:,1,1,1));

%check for spots
if nargin > 1
    SetUserData(imgFigH,spots,0);
    hold on;
    for i=1:size(spots(1).sp,2)
        if strcmp(spots(1).sp(i).type,'spb')
            mark='o';
        elseif strcmp(spots(1).sp(i).type,'kin')
            mark='+';
        else
            mark='r+';
        end;
        if spots(1).sp(i).cord(3)==1
            plot(spots(1).sp(i).cord(1),spots(1).sp(i).cord(2),mark);
        end;
    end;
    if nargin > 2
        histFigH=figure('Name','Histogram');
        SetUserData(imgFigH,mnp,0);
        [hy,hx]=hist(nonzeros(mnp(:,1)),255);
        bar(hx,hy,'stacked');
    end;
end;

%set props
SetUserData(imgFigH,img,0);
stackslideH=findobj(contFigH,'Tag','stackslider');
stackeditH=findobj(contFigH,'Tag','stackedit');
timeslideH=findobj(contFigH,'Tag','timeslider');
timeeditH=findobj(contFigH,'Tag','timeedit');
%set stackslider and edit values
set(stackslideH,'Min',1);
set(stackslideH,'Max',size(img,3));
set(stackslideH,'SliderStep',[1/size(img,3) 10/size(img,3)]);
set(stackslideH,'Value',1);
set(stackeditH,'String','1');

%set timeslider and edit values
set(timeslideH,'Min',0.9999999);
set(timeslideH,'Max',size(img,5));
set(timeslideH,'SliderStep',[1/size(img,5) 10/size(img,5)]);
set(timeslideH,'Value',1);
set(timeeditH,'String','1');

