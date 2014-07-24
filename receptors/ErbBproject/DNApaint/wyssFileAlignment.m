function [handles] = wyssFileAlignment(List,List2,name,varargin)
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
ip.addOptional('thres',10000,@isnumeric);

ip.parse(List,List2,name,varargin{:});

thres = ip.Results.thres;

Imsize = [256,256];
difLim = 2.1075;

%Load data and prepare for processing


for i = 1:numel(List)
    drift = dlmread(List2{i},'\t');   
    save(List{i},'drift','-append');
end

clear drift;

%Image Cell array will contain images for picking drift marker
ImArr = cell([numel(List),1]);

PointList=[];

for i = 1:numel(List);

load(List{i});
%drift correct using trace
pnts=trace(:,2:3)-drift(trace(:,1),:);
pnts = [pnts,trace(:,10)];

ind = pnts(:,3)>thres;

ImArr{i} = hist3(pnts(ind,1:2),'Edges',{[0:0.1:Imsize(1)],[0:0.1:Imsize(2)]});
%ImArr{i} = hist3(pnts(ind,1:2),'Edges',{[0:0.2:Imsize(1)],[0:0.2:Imsize(2)]});
%ImArr{i} = hist3(pnts(:,1:2),'Edges',{[0:Imsize(1)],[0:Imsize(2)]});
PointList = vertcat(PointList, {struct('pnts',pnts(:,1:2),'name',List{i},'fullData',trace,'drift',drift,'dmark',[],'shift',[],'com',[],'ind',ind)});

end

% makes object for figure
f = figure('ToolBar','none');
hp1 = uipanel(f,'Position',[0,0,1,0.8]);
hp2 = uipanel(f,'Position',[0,0.8,1,0.2]);
a = axes('Parent',hp1,'DataAspectRatio',[1,1,1]);
a = get(hp1,'Children');
xslide = uicontrol(hp1,'Style','slider','Min',-100,'Max',100, ...
              'Value',0,'Position',[10,50,30,150], ...
              'Callback',{@shiftadjust}, 'SliderStep',[0.005,0.05]);
yslide = uicontrol(hp1,'Style','slider','Min',-100,'Max',100, ...
              'Value',0,'Position',[10,10,150,30], ...
              'Callback',{@shiftadjust}, 'SliderStep',[0.005,0.05]);
text = uicontrol(hp1,'Style','text','String','x = 0.0  y = 0.0', ...
                 'Position',[10,400,150,30]);
s = uicontrol(hp2,'Style','slider','Min',1,'Max',numel(List), ...
              'Value',1,'Position',[10,10,150,30], ...
              'Callback',{@slideradjust}, 'SliderStep',[1/numel(List),10/numel(List)]);
t = uicontrol(hp2,'Style','text','Position',[10,30,80,20],'String', 'N = 1');
t2 = uicontrol(hp2,'Style','text','Position',[315,10,35,20],'String', 'X lims');
xlim1 = uicontrol(hp2,'Style','edit','Position',[350,10,50,20],'String', '0');
xlim2 = uicontrol(hp2,'Style','edit','Position',[400,10,50,20],'String', '256');
t3 = uicontrol(hp2,'Style','text','Position',[315,30,35,20],'String', 'Y lims');
ylim1 = uicontrol(hp2,'Style','edit','Position',[350,30,50,20],'String', '0');
ylim2 = uicontrol(hp2,'Style','edit','Position',[400,30,50,20],'String', '256');


b = uicontrol(hp2,'Style','checkbox','Position',[10,55,55,20],...
                    'String','Compare','Callback',{@checkboxAction});
push = uicontrol(hp2,'Style','pushbutton','String','GetPnts',...
                    'Position',[180,10,60,20],'Callback',{@selectPnts});
savePL =  uicontrol(hp2,'Style','pushbutton','String','S & Q',...
                    'Position',[180,35,60,20],'Callback',{@savePntList});
check = uicontrol(hp2,'Style','pushbutton','String','Check',...
                    'Position',[250,10,60,45],'Callback',{@checkAlign});

zoom(a,'on');                
h=zoom;
set(h, 'ActionPostCallback', @(x) mypostcallback(h));
                
%set state variables
state = 1; %starting position for display

%set gui data for object usage
setappdata(f,'data',struct('PL',PointList,'name',name));

slider = struct('a',a,'t',t,'b',b,'xs',xslide,'ys',yslide,'i',state,'xlim',[xlim1,xlim2],'ylim',[ylim1,ylim2]);
slider.data = ImArr;
slider.dmark = repmat({[inf,inf]},[size(ImArr,1),1]);
shift = repmat({struct('calc',[0,0],'adj',[0,0],'total',[0,0])},[numel(ImArr),1]);
slider.shift = shift;
setappdata(s,'slider',slider);

shiftad = struct('a',a,'s',s,'xs',xslide,'ys',yslide,'shift',[0,0],'text',text);
setappdata(xslide,'data',shiftad);
setappdata(yslide,'data',shiftad);

pushB = struct('a',a,'s',s);
setappdata(push,'data',pushB);

setappdata(b,'data',struct('s',s));

setappdata(savePL,'data',struct('s',s,'f',f,'name',name));

setappdata(check,'data',struct('s',s,'f',f));

setappdata(h,'data', struct('xlim',[xlim1,xlim2],'ylim',[ylim1,ylim2]));

%creates output structure;
handles = struct('figure', f,'axes',a,'slider',s, 'text',t);

%itialize plot
figure(f);
tmp =slider.data{state};
m = max(log(tmp(:)));
imshow(log(tmp),[0,m]);
colormap('hot');
   
end


function checkAlign(hObject, eventData, handles)
        st = getappdata(hObject,'data');
        slid = getappdata(st.s,'slider');
        d = getappdata(st.f,'data');
        
        i=slid.i;
        
        z=axis;
        
        poly = [z([1,3]);z([2,3]);z([2,4]);z([1,4]);z([1,3])]/10;
                
        
        temp = d(1).PL.pnts(d(1).PL.ind,1:2);
        refpnt = temp(inpolygon(temp(:,1),temp(:,2),poly(:,2),poly(:,1)),:);
        temp = d(i).PL.pnts(d(i).PL.ind,1:2);
        pnt = temp(inpolygon(temp(:,1),temp(:,2),poly(:,2),poly(:,1)),:);
        pnt = pnt - repmat(slid.shift{i}.total,[numel(pnt(:,1)),1]);
        
        scatter(refpnt(:,1),refpnt(:,2),'m.')
        hold;
        scatter(pnt(:,1),pnt(:,2),'g.');
        hold;
        axis(z([3,4,1,2])/10);
end


function slideradjust(hObject, eventData, handles)
    %get data stored in slider
    st = getappdata(hObject,'slider');    
    imgArr = st.data;
    st.i = fix(get(hObject,'Value'));
    axes(st.a);
    shift =st.shift{st.i}.adj*10;
    set(st.xs,'Value',shift(1));
    set(st.ys,'Value',shift(2));
    z = axis;
    set(st.ylim(1),'String',num2str(z(3)/10))
    set(st.ylim(2),'String',num2str(z(4)/10))
    set(st.xlim(1),'String',num2str(z(1)/10))
    set(st.xlim(2),'String',num2str(z(2)/10))
    
    if (get(st.b,'value' ) & ~isnan(st.dmark{1}))
       tImg = repmat(fix(0.5*imgArr{1}),[1,1,3]);
       siz = size(tImg);
       % create a shifted image
       shift = fix((st.shift{st.i}.total)*10);
       ashift = abs(shift);
       temp = zeros(siz(1:2)+2*ashift);
       temp(ashift(1)+1-shift(1):siz(1)+ashift(1)-shift(1),ashift(2)+1-shift(2):siz(2)+ashift(2)-shift(2))= imgArr{st.i};
       tImg(:,:,2)= temp(abs(shift(1))+1:end-abs(shift(1)),abs(shift(2))+1:end-abs(shift(2)));
       
       imshow(tImg,[]);       
       axis(z); 
       colormap('hot');
    else
        tImg = imgArr{st.i};
        m = max(log(tImg(:)));
        imshow(log(tImg),[0,m]);      
        axis(z);
        colormap('hot');
    end
 
    set(st.t,'String',['N = ',num2str(st.i,'%u')]);
    
    %update data stored in slider
    setappdata(hObject,'slider',st);
end

function shiftadjust(hObject,eventData,handles)

st = getappdata(hObject,'data');
x=fix(get(st.xs,'value'));
y=fix(get(st.ys,'value'));
slider = getappdata(st.s,'slider');
i = slider.i;
slider.shift{i}.adj = [x,y]/10;
slider.shift{i}.total = slider.shift{i}.calc+slider.shift{i}.adj;
set(st.text,'String',['x = ',num2str(x/10,'%3.1f'),' y = ',num2str(y/10,'%3.1f')]);
setappdata(st.s,'slider',slider);
slideradjust(st.s, eventData, st.s);

end


function checkboxAction(hObject,eventData,handles)
   st = getappdata(hObject,'data');
   slideradjust(st.s, eventData, st.s)
end

function selectPnts(hObject,eventData,handles)

 data = getappdata(hObject,'data');
 slider = getappdata(data.s,'slider');
 zoom(data.a,'off')
 [x,y] = getpts(data.a);
 zoom(data.a,'on')
 
 
 x=fix(x);
 y=fix(y);
 hold
 width = 51;
 hw = 25;
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
         temp.A =1;
         display(['failed to fit ', num2str(i),'th mark: see red xs']);
         hold; scatter(x(i),y(i),'rx');hold;  
    end
    end
%dmark(i,:)=[(x(i)+(temp.x-(hw+1)))/5,(y(i)+(temp.y-(hw+1)))/5];
dmark(i,:)=[(y(i)+(temp.y-(hw+1)))/10,(x(i)+(temp.x-(hw+1)))/10];
end

dmark(isnan(dmark(:,1)),:)=[];
slider.dmark{slider.i}=dmark;

i = slider.i;

if isfinite(slider.dmark{1}) & ~isempty(slider.dmark{i})
    slider.shift{i}.calc=(slider.dmark{i}-slider.dmark{1});
    slider.shift{i}.total = slider.shift{i}.calc+slider.shift{i}.adj;
else
    display('select drift marker in first image first');
end

setappdata(data.s,'slider',slider);
if (any(isnan(dmark)))
    display('Error in fitting');
end

end


function savePntList(hObject,eventData,handles)

    data = getappdata(hObject,'data');
    slider = getappdata(data.s,'slider');
    fig = getappdata(data.f,'data');

    name = data.name;
    
    PointList = vertcat({fig.PL});
    dmarks = slider.dmark;
    
    for i = 1:numel(PointList)
      PointList{i}.dmark = dmarks{i};
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
%      test = PointList{j}.dmark;
%      [shift,transform,A,B] = driftMarkerRegistration(test,ref,Imsize,difLim);
%      TotalShift = shift - transform.trans;
%      C = test - repmat(TotalShift,[numel(test(:,1)),1]);

     TotalShift = slider.shift{j}.calc+slider.shift{j}.adj;
     PointList{j}.shift = struct('transform',[],'preshift',[],'A',[],'B',[],'Postshift',[],'TotalShift',TotalShift);
     
     %shifts points to align
     tmp = PointList{j}.pnts(:,1:2)-repmat(TotalShift,[numel(PointList{j}.pnts(:,1)),1]);
     PointList{j}.pnts=tmp(~isnan(tmp(:,1)),:);
     
end

 %Since all other Images are shifted to match this one, here only NaNs are
 %removed
 tmp = PointList{j_ref}.pnts(:,1:2);
 PointList{j_ref}.pnts=tmp(~isnan(tmp(:,1)),:);
 
 
save([name,'_manualAlign_PointList.mat'],'PointList');
    
close(data.f)    


end


function mypostcallback(h,eventData)
    st = getappdata(h,'data');
    z = axis;
    z=fix(z/10);
    set(st.ylim(1),'String',num2str(z(3)))
    set(st.ylim(2),'String',num2str(z(4)))
    set(st.xlim(1),'String',num2str(z(1)))
    set(st.xlim(2),'String',num2str(z(2)))
    
end


