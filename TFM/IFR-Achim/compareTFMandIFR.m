%nargin=0;
nargin=0

if nargin < 1 || isempty(TFMforceField)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select the TFM forceField.mat to be used');
       fileStruct=load([pathname filesep filename]);
       TFMforceField=fileStruct.forceField;
end

%read in stack of flow files:
if nargin < 2 || isempty(inputFileList)
   [filename, pathname] = uigetfile({'*.mat';'*.*'}, ...
       'Select the first IFR forceField.mat file to be analyzed');
   
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

n = numel(inputFileList);


for i=1:3
    fileStruct=load(inputFileList{i});
    IFRforceField=fileStruct.forceField;

    figure(i); 
    quiver(IFRforceField.p(:,1),IFRforceField.p(:,2),IFRforceField.f(:,1),IFRforceField.f(:,2),'r')
    hold on;
    quiver(TFMforceField(i).pos(:,1),TFMforceField(i).pos(:,2),TFMforceField(i).vec(:,1),TFMforceField(i).vec(:,2),'b');
    set(gca,'YDir','reverse')
    hold off;  
    
    iTFMforceField(i).pos=IFRforceField.p;
    iTFMforceField(i).vec=zeros(length(IFRforceField.f),2);
    iTFMforceField(i).vec(:,1)=griddata(TFMforceField(i).pos(:,1),TFMforceField(i).pos(:,2),TFMforceField(i).vec(:,1),iTFMforceField(i).pos(:,1),iTFMforceField(i).pos(:,2));
    iTFMforceField(i).vec(:,2)=griddata(TFMforceField(i).pos(:,1),TFMforceField(i).pos(:,2),TFMforceField(i).vec(:,2),iTFMforceField(i).pos(:,1),iTFMforceField(i).pos(:,2));  
    
    figure(100+i); 
    quiver(IFRforceField.p(:,1),IFRforceField.p(:,2),IFRforceField.f(:,1),IFRforceField.f(:,2),'r')
    hold on;
    quiver(iTFMforceField(i).pos(:,1),iTFMforceField(i).pos(:,2),iTFMforceField(i).vec(:,1),iTFMforceField(i).vec(:,2),'b');
    set(gca,'YDir','reverse')
    hold off;
    
    %compare magnitude of the force vectors:
    iTFMforceField(i).mag=sqrt(sum(iTFMforceField(i).vec.*iTFMforceField(i).vec,2));
    %IFRforceField.mag =sqrt(sum( IFRforceField.f  .* IFRforceField.f,2));
    
    figure(1000+i)
    plot(IFRforceField.fMag,iTFMforceField(i).mag,'*');
    
    %compare the angle between the two force vectors:
    
    alphas{i}=180/pi*acos(sum(iTFMforceField(i).vec.*IFRforceField.f,2)./(iTFMforceField(i).mag.*IFRforceField.fMag))
    
    figure(2000+i)
    hist(alphas{i});    
    %alpha{i}=acos(dot(-resForce{i}.Cell{1}.vec,resForce{i}.Cell{2}.vec)/(absForce{i}.Cell{1}*absForce{i}.Cell{2}));
end


