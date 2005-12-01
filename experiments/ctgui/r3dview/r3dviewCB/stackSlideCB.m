function stackSlideCB
% callback for trackViewer

imgFigureH=findobj('Type','figure','Name','View Panel');
if isempty(imgFigureH)
   return;
end;
img=GetUserData(imgFigureH,'img');
spots=GetUserData(imgFigureH,'spots');

val=round(get(gcbo,'Value'));
set(gcbo,'Value',val);
stackeditH=findobj(gcbf,'Style','edit','Tag','stackedit');
timesliderH=findobj(gcbf,'Style','slider','Tag','timeslider');
tp=get(timesliderH,'Value');

set(stackeditH,'String',num2str(val));

%display z-slice
imgH=findobj(imgFigureH,'Type','image');
figure(imgFigureH);
set(imgH,'CData',img(:,:,val,1,tp));
if ~isempty(spots)
    for i=1:size(spots(tp).sp,2)
        if strcmp(spots(tp).sp(i).type,'spb')
            mark='yo';
        elseif strcmp(spots(tp).sp(i).type,'kin')
            mark='y+';
        else
            mark='r+';
        end;
        if round(spots(tp).sp(i).cord(3))==val
            plot(spots(tp).sp(i).cord(1),spots(tp).sp(i).cord(2),mark);
        end;
    end;  
end;

