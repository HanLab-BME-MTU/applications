%this is an example of how to read and use information saved by ConeTrack
%by Denis Tsygankov 

load('save_example.mat');

% detection and tracking parameters:
% critical radius
Ocr=fw;
% neck width
Wcr=FL;
% critical distance
Dcr=fd;
% time gap
Gcr=fg;
% lifetime filter
Tcr=ff;

% time duration of the record
NUMfrms=length(nCx);

% display cone boundaries
fig=figure('Position',[200 200 800 400],'Name','An example of the output reading');
sb=subplot(1,2,1);
hold on;

for fr=1:5:NUMfrms
    Xbnd=nConeBnd{fr}(:,1);
    Ybnd=nConeBnd{fr}(:,2);
    fill3(Xbnd,Ybnd,fr*ones(length(Xbnd),1),'c','FaceAlpha',0.5);
end
% display centroid track
plot3(nCx,nCy,1:NUMfrms,'r','LineWidth',2);
axis equal;
axis ij;
view(15,20);
set(sb,'Box','on');
zlabel('time (frames)');

% display cone's velocity
Vx=nCx(2:end)-nCx(1:(end-1));
Vy=nCy(2:end)-nCy(1:(end-1));
Vel=sqrt(Vx.^2+Vy.^2);
sVel=GFilterA(Vel,20);
sb=subplot(2,2,2);
hold on;
plot(Vel,'b');
plot(sVel,'r','LineWidth',2);
xlabel('time (frames)');
ylabel('velocity (pixel/frame)');
set(sb,'Box','on','XLim',[0 length(nCx)],'XGrid','on','YGrid','on');

% display cone's distance to the terminal point
Dis=sqrt((nCx-nCx(end)).^2+(nCy-nCy(end)).^2);
sb=subplot(2,2,4);
hold on;
plot(Dis,'r','LineWidth',2);
xlabel('time (frames)');
ylabel('distance to the end (pixel)');
set(sb,'Box','on','XLim',[0 length(nCx)+1],'XGrid','on','YGrid','on');



