%Select the region by drawing a polygon where we want to analyze the vector
% field and identify the force distribution.
% Always start with the left-lower corner when the polygon is viewed as a
% deformed rectangle.

figure(gcf); hold off;

load([resultPath 'dispField']);
imshow(rawI{1},[]); hold on;

%Plot the raw data points and displacements and the filtered displacements.
quiver(dataPx{1},dataPy{1},rawDataU1{1}*2,rawDataU2{1}*2,0,'y'); 
quiver(dataPx{1},dataPy{1},sDataU1{1}*15,sDataU2{1}*15,0,'r'); 
plot(fieldPGx,fieldPGy,'go');
plot(fieldPGx,fieldPGy,'b');

%Choose the region for identification.
choice = input(['What is the region for identification?\n' ...
   '  ''0'': keep the old choice.\n' ...
   '  ''1'': the current polygon.\n' ...
   '  ''2'': draw a new polygon.\n' ...
   'Your answer: ']);

if choice == 2
   [bw,recPGx,recPGy] = roipoly;
   plot(recPGx,recPGy,'b');
   plot(recPGx(1),recPGy(1),'go');
elseif choice == 1
   recPGx = fieldPGx;
   recPGy = fieldPGy;
elseif choice == 0
   load([resultPath 'recGeom']);
   plot(recPGx,recPGy,'b');
end

if choice == 1 | choice == 2
   for k = 1:length(recPGx)
      text(recPGx(k),recPGy(k),num2str(k),'color','r');
   end

   %Ask for the indices of the other three vertices of the polygon, 'recPG'  
   % when it is viewed as a deformed rectangle. The first vertex is always 1.
   recPGVI = ones(1,4);
   recPGVI(2:4) = input(['What are the indices of the four vertices ' ...
      'when the polygon is viewed as a deformed rectangle?\n  ' ...
      'Start with the left-lower corner (index: 1) and\n  ' ...
      'enter the other three by going clockwise:']);
   plot(recPGx(recPGVI),recPGy(recPGVI),'go');

   save([resultPath 'recGeom'],'recPGx','recPGy','recPGVI');
end

%Draw the four edges.
curvL = [recPGx(recPGVI(1):recPGVI(2)) recPGy(recPGVI(1):recPGVI(2))].';
curvB = [recPGx(end:-1:recPGVI(4)) recPGy(end:-1:recPGVI(4))].';
curvR = [recPGx(recPGVI(4):-1:recPGVI(3)) ...
recPGy(recPGVI(4):-1:recPGVI(3))].';
curvT = [recPGx(recPGVI(2):recPGVI(3)) recPGy(recPGVI(2):recPGVI(3))].';
plot(curvL(1,:),curvL(2,:),'g');
plot(curvT(1,:),curvT(2,:),'g');
plot(curvR(1,:),curvR(2,:),'g');
plot(curvB(1,:),curvB(2,:),'g');

%Choose the polygon of the body force domain.
choice = input(['What is the domain of the external force?\n' ...
   '  ''-1'': no external force.\n' ...
   '  ''0'': keep the old choice.\n' ...
   '  ''1'': the whole identification region.\n' ...
   '  ''2'': draw a new polygon.\n' ...
   'Your answer: ']);
if choice == -1
   bfDomPGx  = [];
   bfDomPGy  = [];
   bfDomPGVI = [];
elseif choice == 1
   bfDomPGx = recPGx; 
   bfDomPGy = recPGy; 
   bfDomPGVI = recPGVI;
elseif choice == 2
   [bw,bfDomPGx,bfDomPGy] = roipoly;
   plot(bfDomPGx,bfDomPGy,'b');
   for k = 1:length(bfDomPGx)
      text(bfDomPGx(k),bfDomPGy(k),num2str(k),'color','r');
   end
   %plot(bfDomPGx,bfDomPGy,'ro');
   plot(bfDomPGx(1),bfDomPGy(1),'go');

   bfDomPGVI = ones(1,4);
   bfDomPGVI(2:4) = input(['What are the indices of the four vertices ' ...
      'when the polygon is viewed as a deformed rectangle?\n  ' ...
      'Start with the left-lower corner (index: 1) and\n  ' ...
      'enter the other three by going clockwise:']);
   plot(bfDomPGx(bfDomPGVI),bfDomPGy(bfDomPGVI),'go');

   %[bfDomPGx,bfDomPGy] = rotate2D(angle,bfDomPGx-ctrX,-bfDomPGy+ctrY);
end

save([modelPath 'recGeom'],'recPGx','recPGy','recPGVI', ...
   'bfDomPGx','bfDomPGy','bfDomPGVI');
hold off;
