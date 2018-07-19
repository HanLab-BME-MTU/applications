function timeSlideCB
% callback for r3dViewer

imgFigureH=findobj('Type','figure','Name','View Panel');
if isempty(imgFigureH)
   return;
end;
img=GetUserData(imgFigureH,'img');
spots=GetUserData(imgFigureH,'spots');
mnp=GetUserData(imgFigureH,'mnp');

val=round(get(gcbo,'Value'));
set(gcbo,'Value',val);
timeeditH=findobj(gcbf,'Style','edit','Tag','timeedit');
set(timeeditH,'String',num2str(val));

stackliderH=findobj(gcbf,'Style','slider','Tag','stackslider');
sp=get(stackliderH,'Value');

%display z-slice
%imgH=findobj(imgFigureH,'Type','image');
%set(imgH,'CData',img(:,:,sp,1,val));
figure(imgFigureH);
cla;
%uiViewPanelShowImg(img(:,:,sp,1,val), 0, imgFigureH);
imshow(img(:,:,sp,1,val));
if ~isempty(spots)
    for i=1:size(spots(val).sp,2)
        if strcmp(spots(val).sp(i).type,'spb')
            mark='o';
        elseif strcmp(spots(val).sp(i).type,'kin')
            mark='+';
        else
            mark='r+';
        end;
        if spots(val).sp(i).cord(3)==sp
            plot(spots(val).sp(i).cord(1),spots(val).sp(i).cord(2),mark);
        end;
    end;  
end;
if ~isempty(mnp)
    histFigH=findobj(0,'Type','figure','Name','Histogram');
    if ~isempty(histFigH)
        figure(histFigH);
        [hy,hx]=hist(nonzeros(mnp(:,val)),255);
        bar(hx,hy,'stacked');
    end;
end;  