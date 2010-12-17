function [resForce]=cutOutForceField()
% This function is specific for one old dataset in:
% 2010_02_24_Rosa_eCAD

if nargin < 1
   [filename_forceField, pathname_forceField] = uigetfile({'*.mat';'*.*'}, ...
       'Select forceField.mat to be used');
       %the stage drift Transformation:
       fileStruct=load([pathname_forceField filesep filename_forceField]);
       forceField=fileStruct.forceField;
end

%read in Stack of bead images:
if nargin < 2 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.tif';'*.jpg';'*.png';'*.*'}, ...
       'Select First Ecad-Image');
   
   if ~ischar(filename) || ~ischar(pathname)
       return;
   end
   
   inputFileList = getFileStackNames([pathname filesep filename]);
else
    isValid = 1;
    for i = 1:numel(inputFileList)
        isValid = isValid && exist(inputFileList{i}, 'file');
    end
    if ~isValid
        error('Invalid input files.');
    end
end

targetDirCellCellForce=[pathname_forceField,filesep ,'cellCellForces', filesep];
if ~isdir(targetDirCellCellForce)
    mkdir(targetDirCellCellForce)
end


%the points are:
p1=[110   44];
p2=[1265 998];

% determine if points are above or below the line:
m=(p2(2)-p1(2))/(p2(1)-p1(1));

n = numel(inputFileList);
padZeros=floor(log10(2*n))+1;

for i=1:n

    checkVector=(forceField(i).pos(:,2)-p1(2)<m*(forceField(i).pos(:,1)-p1(1)));

    forceFieldOfFrame{i}.Cell{1}.pos=forceField(i).pos(checkVector,:);
    forceFieldOfFrame{i}.Cell{1}.vec=forceField(i).vec(checkVector,:);

    forceFieldOfFrame{i}.Cell{2}.pos=forceField(i).pos(~checkVector,:);
    forceFieldOfFrame{i}.Cell{2}.vec=forceField(i).vec(~checkVector,:);

    ecad(i).image=double(imread(inputFileList{i}));
    figure(1)
    imagesc(ecad(i).image)
    colormap('gray')
    hold on
    quiver(forceFieldOfFrame{i}.Cell{1}.pos(:,1),forceFieldOfFrame{i}.Cell{1}.pos(:,2),forceFieldOfFrame{i}.Cell{1}.vec(:,1),forceFieldOfFrame{i}.Cell{1}.vec(:,2),'b');
    quiver(forceFieldOfFrame{i}.Cell{2}.pos(:,1),forceFieldOfFrame{i}.Cell{2}.pos(:,2),forceFieldOfFrame{i}.Cell{2}.vec(:,1),forceFieldOfFrame{i}.Cell{2}.vec(:,2),'g');
    set(gca,'YDir','reverse')
    hold off

    % Further constraints which forces to pick:
    % E.g. define an ellipse by two points and the pin-length:
    point1=[562 768];
    point2=[908 327];
    pinL=sqrt(sum((point1-point2).^2))+2*50;

    cell1_point1=repmat(point1,length(forceFieldOfFrame{i}.Cell{1}.pos),1);
    cell1_point2=repmat(point2,length(forceFieldOfFrame{i}.Cell{1}.pos),1);
    cell1_inVec=(sqrt(sum((forceFieldOfFrame{i}.Cell{1}.pos-cell1_point1).^2,2))+sqrt(sum((forceFieldOfFrame{i}.Cell{1}.pos-cell1_point2).^2,2))<pinL);

    cell2_point1=repmat(point1,length(forceFieldOfFrame{i}.Cell{2}.pos),1);
    cell2_point2=repmat(point2,length(forceFieldOfFrame{i}.Cell{2}.pos),1);
    cell2_inVec=(sqrt(sum((forceFieldOfFrame{i}.Cell{2}.pos-cell2_point1).^2,2))+sqrt(sum((forceFieldOfFrame{i}.Cell{2}.pos-cell2_point2).^2,2))<pinL);

    forceFieldConstraint{i}.Cell{1}.pos=forceFieldOfFrame{i}.Cell{1}.pos(cell1_inVec,:);
    forceFieldConstraint{i}.Cell{1}.vec=forceFieldOfFrame{i}.Cell{1}.vec(cell1_inVec,:);
    
    forceFieldConstraint{i}.Cell{2}.pos=forceFieldOfFrame{i}.Cell{2}.pos(cell2_inVec,:);
    forceFieldConstraint{i}.Cell{2}.vec=forceFieldOfFrame{i}.Cell{2}.vec(cell2_inVec,:);
    
    % Center Point for plotting
    center_pt=[744 560];    
    
    % Sum up the force vectors for each cell:
    resForce{i}.Cell{1}.pos=center_pt;
    resForce{i}.Cell{1}.vec=sum(forceFieldConstraint{i}.Cell{1}.vec,1);
    
    resForce{i}.Cell{2}.pos=center_pt;
    resForce{i}.Cell{2}.vec=sum(forceFieldConstraint{i}.Cell{2}.vec,1);

    % Replot the results:
    figure(2)
    imagesc(ecad(i).image)
    colormap('gray')
    hold on
    quiver(forceFieldConstraint{i}.Cell{1}.pos(:,1),forceFieldConstraint{i}.Cell{1}.pos(:,2),forceFieldConstraint{i}.Cell{1}.vec(:,1),forceFieldConstraint{i}.Cell{1}.vec(:,2),'b');
    quiver(resForce{i}.Cell{1}.pos(1),resForce{i}.Cell{1}.pos(2),resForce{i}.Cell{1}.vec(1),resForce{i}.Cell{1}.vec(2),10^(-2),'b');
    
    quiver(forceFieldConstraint{i}.Cell{2}.pos(:,1),forceFieldConstraint{i}.Cell{2}.pos(:,2),forceFieldConstraint{i}.Cell{2}.vec(:,1),forceFieldConstraint{i}.Cell{2}.vec(:,2),'r');
    quiver(resForce{i}.Cell{2}.pos(1),resForce{i}.Cell{2}.pos(2),resForce{i}.Cell{2}.vec(1),resForce{i}.Cell{2}.vec(2),10^(-2),'r');
    set(gca,'YDir','reverse')
    saveas(gcf,[targetDirCellCellForce,'cellCellForce',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
    hold off

    % give rel. error (the !sum! not the difference has to be zero):

    relError{i}.inXYcomp=(resForce{i}.Cell{1}.vec+resForce{i}.Cell{2}.vec)./resForce{i}.Cell{1}.vec;
    absError{i}=(sqrt(sum((resForce{i}.Cell{1}.vec).^2))-sqrt(sum((resForce{i}.Cell{2}.vec).^2)))/sqrt(sum((resForce{i}.Cell{1}.vec).^2));
    absForce{i}.Cell{1}=sqrt(sum((resForce{i}.Cell{1}.vec).^2));
    absForce{i}.Cell{2}=sqrt(sum((resForce{i}.Cell{2}.vec).^2));
    alpha{i}=acos(dot(-resForce{i}.Cell{1}.vec,resForce{i}.Cell{2}.vec)/(absForce{i}.Cell{1}*absForce{i}.Cell{2}));
end

% Solve for the displacement field of the two forcefields and then
% calculate the elastic energy:
for i=1:n
    for j=1:2
        % Start with the first cell:
        %extend the functions to a regular grid:
        
%!!!    %maybe instead of griddata I should use TriScatterData, might be
        %faster:
        force_x=@(x,y) myGriddata(forceFieldConstraint{i}.Cell{j}.pos(:,1),forceFieldConstraint{i}.Cell{j}.pos(:,2),forceFieldConstraint{i}.Cell{j}.vec(:,1),x,y,'cubic');
        force_y=@(x,y) myGriddata(forceFieldConstraint{i}.Cell{j}.pos(:,1),forceFieldConstraint{i}.Cell{j}.pos(:,2),forceFieldConstraint{i}.Cell{j}.vec(:,2),x,y,'cubic');

        % Here we need to calculte the displacements only at the positions where
        % the force field doesn't vanish. For the elastic energy we want to 
        % calculate int[u(x)*f(x)] and the integrand is zero if f=0.
        E=8000;
        
        xmin=min(min(forceFieldConstraint{i}.Cell{j}.pos(:,1)));
        xmax=max(max(forceFieldConstraint{i}.Cell{j}.pos(:,1)));
        ymin=min(min(forceFieldConstraint{i}.Cell{j}.pos(:,2)));
        ymax=max(max(forceFieldConstraint{i}.Cell{j}.pos(:,2)));
        
        [ux uy]=fwdSolution(forceFieldConstraint{i}.Cell{j}.pos(:,1),forceFieldConstraint{i}.Cell{j}.pos(:,2),E,xmin,xmax,ymin,ymax,force_x,force_y,'conv');
        displFieldConstraint{i}.Cell{j}.pos=forceFieldConstraint{i}.Cell{j}.pos;
        displFieldConstraint{i}.Cell{j}.vec(:,1)=ux;
        displFieldConstraint{i}.Cell{j}.vec(:,2)=uy;
        
        factor=forceField(i).par.gridSpacing^2*forceField(i).par.pixSize_mu^3/10^6;        
        elEnergy(i).Cell{j}= 1/2*sum(sum(displFieldConstraint{i}.Cell{j}.vec.*forceFieldConstraint{i}.Cell{j}.vec))*factor;
        ratioCCForceOverElEnergy(i).Cell{j}=absForce{i}.Cell{j}/elEnergy(i).Cell{j};
        clear('ux','uy')
    end
end

save([pathname_forceField, filesep, 'cellCellForces.mat'], 'forceFieldConstraint', 'displFieldConstraint', 'resForce', 'relError', 'absError', 'absForce', 'alpha', 'elEnergy', 'ratioCCForceOverElEnergy');

startPt=1;
endPt=n;

figure(3)
for i=startPt:endPt
    plot(i,relError{i}.inXYcomp(1),'or')
    hold on
    plot(i,relError{i}.inXYcomp(2),'og')
    title('Development of rel. error over time')
end
hold off


figure(4)
for i=startPt:endPt
    plot(i,absError{i},'ob')
    hold on
    title('Development of abs. error over time')
end
hold off

figure(5)
for i=startPt:endPt
    plot(i,alpha{i}*360/(2*pi),'ob')
    hold on
    title('Development of angular deviations over time')
    ylim([0 180])
end
hold off

maximumForce=0;
figure(6)
for i=startPt:endPt
    plot(absForce{i}.Cell{1},absForce{i}.Cell{2},'ob')
    hold on
    title('residual force of cell2 over residual force of cell1  ')
    currentMax=max(absForce{i}.Cell{:});
    if currentMax>maximumForce
        maximumForce=currentMax;
    end
end
plot([0 maximumForce],[0 maximumForce],'--k')
hold off

figure(7)
for i=startPt:endPt
    colorStr='ob';
    for j=1:2
        if j==2
            colorStr='or';
        end
        plot(i,elEnergy(i).Cell{j},colorStr)
        hold on
        title('Elastic energy invested by each cell') 
    end
end
hold off

figure(8)
for i=startPt:endPt
    colorStr='ob';
    for j=1:2
        if j==2
            colorStr='or';
        end
        plot(i,ratioCCForceOverElEnergy(i).Cell{j},colorStr)
        hold on
        title('Ratio of cell-cell-interaction force and elastic energy') 
    end
end
hold off


function [Zout]=myGriddata(Xin,Yin,Zin,Xout,Yout,method)
    Zout=griddata(Xin, Yin, Zin, Xout, Yout, method);
    nanMat=isnan(Zout);
    Zout(nanMat)=0;
end

end


