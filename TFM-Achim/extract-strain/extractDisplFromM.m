%go to the project folder and execute this function:
function [strain]=extractDisplFromM(doCorrectStageDrift,tol)

%doCorrectStageDrift=1;
%tol=1;

%readin the magic position matrix:
load([pwd,'/tack/mpm.mat'])
[rows,cols,frames]=size(M);
%Probably MPM is integer valued if tracking method was chosen not to be
%subpixel resolution.

%refPic = imread([pwd,'/data/bead001.tif']);

k=1;
j=1;
while k<frames+1
    myM=M(:,:,k);
    %remove zeros:
    myM(~isfinite(1./myM(:,1)),:)=[];
    myM(~isfinite(1./myM(:,3)),:)=[];
    
    %the first column in M is the y1 component of the position vector.
    %the second column in M is the x1 component of the position vector.
    %the third column in M is the y2 component of the position vector.
    %the fourth column in M is the x2 component of the position vector.
    strain(j).pos(:,1)=myM(:,2);
    strain(j).pos(:,2)=myM(:,1);
    
    strain(j).vec(:,1)=myM(:,4)-myM(:,2);
    strain(j).vec(:,2)=myM(:,3)-myM(:,1);

    %if j==1 || j==18
    %    figure(j)
    %    imagesc(refPic)
    %    hold on
    %    quiver(strain(j).pos(:,1),strain(j).pos(:,2),strain(j).vec(:,1),strain(j).vec(:,2))
    %    colormap(gray)
    %    hold off
    %end
    
    k=k+1;
    jmax=j;
    j=j+1;
end

% Correct for stage drift using the mpm matrix, stored in /track/mpm.mat.
% This is very basic and should be imporved using the registration method
% implemented by Sylvain.
correctionVector=[];
if doCorrectStageDrift==1
    %define the region where the image should be stationary.
    [~, cols_refPic]= size(refPic);
    xmin=1;
    ymin=1;
    xmax=cols_refPic;
    ymax=150;
    
    corr_vec_x=zeros(length(strain),1);
    corr_vec_y=zeros(length(strain),1);
    
    for j=1:length(strain)
        tempStrain(:,1:2)=strain(j).pos;
        tempStrain(:,3:4)=strain(j).vec;
        
        %delete all values that are not in the image area defined by
        %xmin, xmax, ymin, ymax:
        tempStrain(((tempStrain(:,1)<xmin)+(tempStrain(:,1)>xmax)+(tempStrain(:,2)<ymin)+(tempStrain(:,2)>ymax))>0,:)=[];
        
        if j==1
            figure(100)
            imagesc(refPic)
            hold on
            plot(tempStrain(:,1),tempStrain(:,2),'or')
            hold off
        end
        
        corr_vec_x(j)=-mean(tempStrain(:,3));
        corr_vec_y(j)=-mean(tempStrain(:,4));
        
        if sqrt(corr_vec_x(j).^2+corr_vec_y(j).^2)>tol
            correctedStrain(j).vec(:,1)=strain(j).vec(:,1)+corr_vec_x(j);
            correctedStrain(j).vec(:,2)=strain(j).vec(:,2)+corr_vec_y(j);
            correctedStrain(j).pos=strain(j).pos;
        else
            correctedStrain(j).vec=strain(j).vec;
            correctedStrain(j).pos=strain(j).pos;
        end
        
        clear('tempStrain');
    
    end
    
    %correct only if stage shift is larger than pixel:
    correctionVector=[(sqrt(corr_vec_x.^2+corr_vec_y.^2)>tol).*corr_vec_x (sqrt(corr_vec_x.^2+corr_vec_y.^2)>tol).*corr_vec_y];
end

if sum(sum(abs(correctionVector)))>0
    %for j=1:length(correctedStrain)
    %    figure(10*j)
    %    imagesc(refPic)
    %    hold on
    %    quiver(correctedStrain(j).pos(:,1),correctedStrain(j).pos(:,2),correctedStrain(j).vec(:,1),correctedStrain(j).vec(:,2))
    %    colormap(gray)
    %    hold off
    %end
    'used corrected displacements'
    clear('strain')
    strain=correctedStrain;
else
    'used un-corrected displacements'
end
quiver(strain(1).pos(:,1),strain(1).pos(:,2),strain(1).vec(:,1),strain(1).vec(:,2),0)
save([pwd,'/mech/strain.mat'],'strain')
