%Select the polygon that defines the lamollipodium.
% Always start with the left-lower corner when the polygon is viewed as a
% deformed rectangle.

%Get the image of the cell.
for jj = 1:length(imgFile)
   rawI{jj} = imread(imgFile{jj});
end

%Show the image 
figure(gcf); hold off;
imshow(rawI{1},[]); axis on; hold on;

%Draw the boundary of the field.
choice = input(['Draw the boundary of the field:\n' ...
   '  ''0'': keep the old choice.\n' ...
   '  ''1'': draw a new polygon.\n' ...
   'Your answer: ']);

if choice == 1
   [bw,fieldPGx,fieldPGy] = roipoly;
   save([modelPath 'fieldGeom'], 'fieldPGx', 'fieldPGy');
elseif choice == 0
   load([modelPath 'fieldGeom']);
end

for k = 1:length(fieldPGx)
   text(fieldPGx(k),fieldPGy(k),num2str(k),'color','r');
end
plot(fieldPGx,fieldPGy,'go');
plot(fieldPGx,fieldPGy,'b');

