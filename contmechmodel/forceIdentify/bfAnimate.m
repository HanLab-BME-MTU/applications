%This script file shows the animation of identified body force over time.

%Load the identified body force.
load([resultPath 'bfId']);

%Get the image of the cell.
cellImg = {};
for jj = 1:numTimeSteps
   cellImg{jj} = imread(imgFile{jj});
end

%M = vectorFieldAnimate([bfDisplayPx bfDisplayPy],recBF,100000, ...
%   'bgImg',cellImg,'vColor','r','colorMap','default');
M = vectorFieldAnimate([bfDisplayPx bfDisplayPy],recBF,bfScale, ...
   'bgImg',cellImg,'aviMovie',{[resultPath 'bfMovie'] {'fps' 10} }, ...
   'colorMap','autumn');
