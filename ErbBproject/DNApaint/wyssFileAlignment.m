function [handles] = wyssFileAlignment(List,List2,name)
%function takes in a PointList recieved from the wyss and creates a GUI
%interface for selecting the drift marker for alignment
%
% inputs: 
%        List, A cell array containing the location of .mat files from 
%              the Wyss to be aligned
%        List2, A cell array contain the location of the .txt files 
%               contianing information for drift corrections (must be in an
%               order such that the files are paired to List)
%        name, the name for saving the files
%
%
% Outputs:
%        PointList, (saved) returns the PointList with the added information about
%        the location of the driftMarkers and save it as name.mat
%
%        handles, structure containing the object handles for the gui
%
%
% Always sets the first element of the PointList to be the reference image
% when preview mode is set Image 1 is displayed in purple and the compared
% image is in green
%
%Written 3/17/2014
%by Jeffrey Werbin
%Harvard Medical School
%

ip=inputParser;
ip.CaseSensitive=false;
ip.StructExpand=true;

ip.addRequired('List',@iscell);
ip.addRequired('List2',@iscell);
ip.addRequired('name',@ischar);

Imsize = [256,256];
difLim = 2.1075;

%Load data and prepare for processing


for i = numel(List)
    drift = dlmread(List2{i},'\t');   
    save(List{i},'drift','-append');
end

%Image Cell array will contain images for picking drift marker
ImArr = cell([numel(List),1]);

for i = 1:numel(List);

load(List{i});
%drift correct using trace
pnts=trace(:,2:3)-drift(trace(:,1),:);
pnts = [pnts,trace(:,10)];

PointList=[];
ImArr{i} = hist3(pnts(:,1:2),'Edges',{[0:0.2:Imsize(1)],[0:0.2:Imsize(2)]});
%ImArr{i} = hist3(pnts(:,1:2),'Edges',{[0:Imsize(1)],[0:Imsize(2)]});
PointList = vertcat(PointList, {struct('pnts',pnts(:,1:2),'name',List{i},'fullData',trace,'drift',drift,'dmark',[],'shift',[],'com',[])});

end

% makes object for figure
f = figure('ToolBar','none');
hp1 = uipanel(f,'Position',[0,0,1,0.8]);
hp2 = uipanel(f,'Position',[0,0.8,1,0.2]);
a = axes('Parent',hp1,'DataAspectRatio',[1,1,1]);
a = get(hp1,'Children');
s = uicontrol(hp2,'Style','slider','Min',1,'Max',numel(List), ...
              'Value',1,'Position',[10,10,150,30], ...
              'Callback',{@slideradjust}, 'SliderStep',[1/numel(List),10/numel(List)]);
t = uicontrol(hp2,'Style','text','Position',[10,30,80,20],'String', 'N = 1');
b = uicontrol(hp2,'Style','checkbox','Position',[10,55,55,20],...
                    'String','Compare','Callback',{@checkboxAction});
push = uicontrol(hp2,'Style','pushbutton','String','GetPnts',...
                    'Position',[180,10,60,20],'Callback',{@selectPnts});
savePL =  uicontrol(hp2,'Style','pushbutton','String','S & Q',...
                    'Position',[180,35,60,20],'Callback',{@savePntList});

%set state variables
state = 1; %starting position for display

%set gui data for object usage
setappdata(f,'data',struct('PL',PointList,'name',name));

slider = struct('a',a,'t',t,'b',b,'i',state);
slider.data = ImArr;
slider.dmark = cell(size(ImArr));
setappdata(s,'slider',slider);

pushB = struct('a',a,'s',s);
setappdata(push,'data',pushB);

setappdata(b,'data',struct('s',s));

setappdata(savePL,'data',struct('s',s,'f',f));

%creates output structure;
handles = struct('figure', f,'axes',a,'slider',s, 'text',t);

%itialize plot
figure(f);
tmp =slider.data{state};
m = max(tmp(:));
imshow(tmp,[0,m]);
   
end





function slideradjust(hObject, eventData, handles)
    %get data stored in slider
    st = getappdata(hObject,'slider');    
    imgArr = st.data;
    st.i = fix(get(hObject,'Value'));
    axes(st.a);
    
    if get(st.b,'value' )
       tImg = repmat(fix(0.5*imgArr{1}),[1,1,3]);
       siz = size(tImg);
       % create a shifted image
       shift = fix((st.dmark{st.i}-st.dmark{1})*5);
       ashift = abs(shift);
       temp = zeros(siz(1:2)+2*ashift);
       temp(ashift(1)+1-shift(1):siz(1)+ashift(1)-shift(1),ashift(2)+1-shift(2):siz(2)+ashift(2)-shift(2))= imgArr{st.i};
       tImg(:,:,2)= temp(abs(shift(1))+1:end-abs(shift(1)),abs(shift(2))+1:end-abs(shift(2)));
       
       imshow(tImg,[]);       
        
    else
        tImg = imgArr{st.i};
        m = max(log(tImg(:)));
        imshow(log(tImg),[0,m]);      
 
    end
 
    set(st.t,'String',['N = ',num2str(st.i,'%u')]);
    
    %update data stored in slider
    setappdata(hObject,'slider',st);
end

function checkboxAction(hObject,eventData,handles)
   st = getappdata(hObject,'data');
   slideradjust(st.s, eventData, st.s)
end

function selectPnts(hObject,eventData,handles)

 data = getappdata(hObject,'data');
 slider = getappdata(data.s,'slider');
 [x,y] = getpts(data.a);
 
 x=fix(x);
 y=fix(y);
 hold
 width = 31;
 hw = 15;
 plot([x-hw,x+hw,x+hw,x-hw,x-hw]',[y-hw,y-hw,y+hw,y+hw,y-hw]','w');
 hold;
 
 
 img = slider.data{slider.i};
 m = max(img(:));
 %There is always at least one
   n = numel(x);
dmark = zeros([n,2]);
for i=1:n
	sigma =1;
    temp = struct('A',NaN);
    while isnan(temp.A)
    temp = fitGaussians2D(img(y(i)-hw:y(i)+hw,x(i)-hw:x(i)+hw),hw+1,hw+1,m,sigma,0,'xyAcs');
    %a bit hacky but if sigma = 1 fails it tries up to sigma = 7 until a
    %fit is reach
    sigma = sigma + 1;
    if sigma >10
        temp.A =1
    end
    end
%dmark(i,:)=[(x(i)+(temp.x-(hw+1)))/5,(y(i)+(temp.y-(hw+1)))/5];
dmark(i,:)=[(y(i)+(temp.y-(hw+1)))/5,(x(i)+(temp.x-(hw+1)))/5];
end

slider.dmark{slider.i}=dmark;

setappdata(data.s,'slider',slider);
if (any(isnan(dmark)))
    display('Error in fitting');
end

end


function savePntList(hObject,eventData,handles)

    data = getappdata(hObject,'data');
    slider = getappdata(data.s,'slider');
    fig = getappdata(data.f,'data');
    
    PointList = fig.PL;
    dmarks = slider.dmark;
    
    for i = numel(PointList)
        PointList(i).dmark = dmarks{i};
    end
    
%aligns images a la exchangePaintAlignment (copied from). Only
%translational alignment

k = 1:numel(PointList);
j_ref = 1;

k(j_ref)=[];

ref = PointList{j_ref}.dmark;
s_ref = size(ref);
PointList{j_ref}.shift = struct('transform',[],'preshift',[],'A',[],'B',[],'Postshift',ref);

for j = k
     test = PointList{j}.dmark;
     [shift,transform,A,B] = driftMarkerRegistration(test,ref,Imsize,difLim);
     TotalShift = shift - transform.trans;
     C = test - repmat(TotalShift,[numel(test(:,1)),1]);
     PointList{j}.shift = struct('transform',transform,'preshift',shift,'A',A,'B',B,'Postshift',C,'TotalShift',TotalShift);
     
     %shifts points to align
     tmp = PointList{j}.pnts(:,1:2)-repmat(TotalShift,[numel(PointList{j}.pnts(:,1)),1]);
     PointList{j}.pnts=tmp(~isnan(tmp(:,1)),:);
     
end

 %Since all other Images are shifted to match this one, here only NaNs are
 %removed
 tmp = PointList{j_ref}.pnts(:,1:2);
 PointList{j_ref}.pnts=tmp(~isnan(tmp(:,1)),:);
 
 
save([name,'.mat'],'PointList');
    
close(data.f)    


end