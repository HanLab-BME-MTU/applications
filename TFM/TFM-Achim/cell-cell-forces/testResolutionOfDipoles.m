function [shift]=testResolutionOfDipoles(displField)
meshPtsFwdSol=2^11;
pixSize_mu=0.163;

% to test this function run this code:
x0=displField(1).pos(:,1);
y0=displField(1).pos(:,2);

[xPlot,yPlot]=meshgrid(1:1:1000,1:1:1000);

origin_init(1)=(max(x0(:))-min(x0(:)))/2;
origin_init(2)=(max(y0(:))-min(y0(:)))/2;

s=1000

% in pixel:
r=round(1/pixSize_mu)

% in pixel:
d=round(30/pixSize_mu)

toDoListSep=0:0.5:5;
toDoListOrigin=rand(50,2)*100;
% in pixel:
for j=1:length(toDoListOrigin)
    for k=1:length(toDoListSep)
        muVal=toDoListSep(k);

        sep=round(muVal/pixSize_mu);
        origin(1)=origin_init(1)+toDoListOrigin(j,1);
        origin(2)=origin_init(2)+toDoListOrigin(j,2);
        origin
        display(['Assumed 0-spacing between neighboring adhesions in mu= ',num2str(sep*pixSize_mu)])


        [fx fy]=twoOppFdipoles(xPlot,yPlot,origin,s,r,d,sep);
%         figure(1)
%         title('Assumed force field')
%         quiver(xPlot,yPlot,fx,fy)

        yModu_Pa=10^4;
        force_x=@(x,y) twoOppFdipoles(x,y,origin,s,r,d,sep);
        force_y=@(x,y) zeros(size(x));

        [ux uy x_grid y_grid]=fwdSolution(x0,y0,yModu_Pa,[],[],[],[],force_x,force_y,'fft',[],meshPtsFwdSol);

%         figure(2)
%         title('Calculated displacement Field')
%         quiver(x_grid,y_grid,ux,uy)

        % here begins the reconstruction of the force Field:
        displFieldModel(1).pos(:,1)=x_grid(:);
        displFieldModel(1).pos(:,2)=y_grid(:);
        displFieldModel(1).vec(:,1)=ux(:);
        displFieldModel(1).vec(:,2)=uy(:);

        regParam = 10^(-6);
        method='FTTC';
        pRatio=0.5;

        [reg_grid,~,~,gridSpacing]=createRegGridFromDisplField(displFieldModel);

        for i=1:length(displFieldModel)
            [grid_mat,iu_mat, i_max,j_max] = interp_vec2grid(displFieldModel(i).pos, displFieldModel(i).vec,[],reg_grid);
            %[grid_mat,u, i_max,j_max] = interp_vec2grid(myStrain(1).pos, myStrain(1).vec,gridSpacing, [])

            %figure(3) 
            %surf(grid_mat(:,:,1),grid_mat(:,:,2),u(:,:,1));

            %figure(4) 
            %surf(grid_mat(:,:,1),grid_mat(:,:,2),u(:,:,2));
            % Now calculte the results:

            if strcmp(method,'FastBEM')
                % If grid_mat=[], then an optimal hexagonal force mesh is created
                % given the bead locations defined in displFieldModel:
                tic
                [pos_f, force, forceMesh, M, pos_u, u, sol_coef]=reg_FastBEM_TFM(grid_mat, displFieldModel, i, yModu_Pa, pRatio, regParam, meshPtsFwdSol);
                display('The total time for calculating the FastBEM solution: ')
                toc

                % The following values are only stored for the BEM-method.
                forceField(i).par.forceMesh = forceMesh;
                forceField(i).par.sol_coef  = sol_coef;
                forceField(i).par.M         = M;
                forceField(i).par.pos       = pos_u;
                forceField(i).par.u         = u;
                forceField(i).par.meshPtsFwdSol=meshPtsFwdSol;
            else
                [pos_f,~,force,~,~,~] = reg_fourier_TFM(grid_mat, iu_mat, yModu_Pa, pRatio, pixSize_mu, gridSpacing, i_max, j_max, regParam);
            end   

%             figure(100)
%             quiver(pos_f(:,1),pos_f(:,2),force(:,1),force(:,2))
%         %    xlim([1 theXlim])
%         %    ylim([1 theYlim])
%         %    set(gca,'YDir','reverse')%,'XTick',[],'YTick',[])
%             title(['Force field frame no: ',num2str(i)])
%             %saveas(gcf,[targetDirForce,'forceField',num2str(i,['%0.',int2str(padZeros),'d']),'.tiff'],'tiffn');
%             set(gca, 'DataAspectRatio', [1,1,50],'YDir','reverse')%,'XTick',[],'YTick',[])


            % Fill in the values to be stored:
            forceField(i).pos=pos_f;
            forceField(i).vec=force;
            forceField(i).posShifted=[]; % this will be calculated when needed
            forceField(i).vecReIntp =[]; % this will be calculated when needed
            forceField(i).par.yModu_Pa   =yModu_Pa;
            forceField(i).par.pRatio     =pRatio;
            forceField(i).par.pixSize_mu =pixSize_mu;
            forceField(i).par.regParam   =regParam;
            forceField(i).par.gridSpacing=gridSpacing;
            forceField(i).par.method     =method;
            forceField(i).doc='Correct pairings: (pos,vec):raw, (posShifted,vec):shifted according to displ, (pos,vecReIntp):reinterp to original grid';
        end

        % force on each adhesion:
        A=pi*r^2/gridSpacing^2;
        f_GT=A*s;

        % boundary offset:
        offset=30;

        % sum all the forces left and right of the origin:
        checkLeft_x=((forceField(1).pos(:,1)<origin(1)) + (forceField(1).pos(:,1)>(origin(1)-(sep/2+2*r+d+offset))))==2;
        checkLeft_y=((forceField(1).pos(:,2)<origin(2)+offset) + (forceField(1).pos(:,2)>(origin(2)-offset)))==2;
        checkLeft_xny=(checkLeft_y + checkLeft_x)==2;
        forceFieldLeft.pos=forceField(1).pos(checkLeft_xny,:);
        forceFieldLeft.vec=forceField(1).vec(checkLeft_xny,:);

        checkRight_x=((forceField(1).pos(:,1)>=origin(1)) + (forceField(1).pos(:,1)<(origin(1)+(sep/2+2*r+d+offset))))==2;
        checkRight_y=((forceField(1).pos(:,2)<origin(2)+offset) + (forceField(1).pos(:,2)>(origin(2)-offset)))==2;
        checkRight_xny=(checkRight_y + checkRight_x)==2;
        forceFieldRight.pos=forceField(1).pos(checkRight_xny,:);
        forceFieldRight.vec=forceField(1).vec(checkRight_xny,:);

        shift(j).f_left(k,:) =sum( forceFieldLeft.vec,1);
        shift(j).f_right(k,:)=sum(forceFieldRight.vec,1);

        % relative error to the measured Traction field:
        shift(j).f_left_rel_err(k,:)=shift(j).f_left(k,:)./sum( abs(forceFieldLeft.vec),1);
        shift(j).f_right_rel_err(k,:)=shift(j).f_right(k,:)./sum( abs(forceFieldRight.vec),1);

        % relative to ground truth magnitude:
        shift(j).f_left_rel_err_GT(k,:)=shift(j).f_left(k,:)./[2*f_GT,0];
        shift(j).f_right_rel_err_GT(k,:)=shift(j).f_right(k,:)./[2*f_GT,0];

        % relative error of the sum of magnitudes:
        shift(j).sum_abs_val_left_rel_error(k,:)=sum(abs(forceFieldLeft.vec),1)/[2*f_GT,0];
        shift(j).sum_abs_val_right_rel_error(k,:)=sum(abs(forceFieldRight.vec),1)/[2*f_GT,0];

%         figure(101)
%         quiver(forceFieldLeft.pos(:,1),forceFieldLeft.pos(:,2),forceFieldLeft.vec(:,1),forceFieldLeft.vec(:,2),'r')
%         hold on;
%         quiver(forceFieldRight.pos(:,1),forceFieldRight.pos(:,2),forceFieldRight.vec(:,1),forceFieldRight.vec(:,2),'g')
%         [fx fy]=twoOppFdipoles(forceField(1).pos(:,1),forceField(1).pos(:,2),origin,s,r,d,sep);
%         quiver(forceField(1).pos(:,1),forceField(1).pos(:,2),fx,fy,'k')
%         hold off;

        % Parameters:
        shift(j).meshPtsFwdSol=meshPtsFwdSol;
    end
end

sep_um=round(toDoListSep/pixSize_mu)*pixSize_mu;
combined_x=horzcat(shift.sum_abs_val_left_rel_error);
combined_x(:,2:2:end)=[];
ave_sum_abs_val_left_rel_error=mean(abs(combined_x),2);
std_sum_abs_val_left_rel_error=std(abs(combined_x),[],2);
display('seperation [um], rel. error., std')
horzcat(sep_um',ave_sum_abs_val_left_rel_error,std_sum_abs_val_left_rel_error)


combined_x=horzcat(shift.f_left_rel_err_GT);
combined_x(:,2:2:end)=[];
ave_f_left_rel_err_GT=mean(abs(combined_x),2);
std_f_left_rel_err_GT=std(abs(combined_x),[],2);
display('seperation [um], rel. error., std')
horzcat(sep_um',ave_f_left_rel_err_GT,std_f_left_rel_err_GT)

combined_x=horzcat(shift.f_left_rel_err);
combined_x(:,2:2:end)=[];
ave_f_left_rel_err=mean(abs(combined_x),2);
std_f_left_rel_err=std(abs(combined_x),[],2);
display('seperation [um], rel. error., std')
horzcat(sep_um',ave_f_left_rel_err,std_f_left_rel_err)


% s=1000; r=round(1/pixSize_mu); d=round(30/pixSize_mu); regParam =
% 10^(-6); offset=30;
% Yields:
% seperation [um], rel. error., std
%          0    0.3073    0.0609
%     0.4890    0.2352    0.0552
%     0.9780    0.1821    0.0501
%     1.4670    0.1425    0.0472
%     1.9560    0.1114    0.0458
%     2.4450    0.0856    0.0427
%     2.9340    0.0688    0.0369
%     3.4230    0.0559    0.0345
%     4.0750    0.0393    0.0311
%     4.5640    0.0304    0.0280
%     5.0530    0.0257    0.0252


% s=1000; r=round(1/pixSize_mu); d=round(30/pixSize_mu); regParam =
% 10^(-7); offset=30;
% Yields:
% seperation [um], rel. error., std
%          0    0.1317    0.0991
%     0.4890    0.0918    0.0800
%     0.9780    0.0742    0.0616
%     1.4670    0.0579    0.0541
%     1.9560    0.0472    0.0467
%     2.4450    0.0399    0.0405
%     2.9340    0.0327    0.0357
%     3.4230    0.0265    0.0297
%     4.0750    0.0218    0.0224
%     4.5640    0.0184    0.0170
%     5.0530    0.0178    0.0148


% s=1000; r=round(1/pixSize_mu); d=round(30/pixSize_mu); regParam =
% 10^(-8); offset=30;
% Yields: (Similar values for ground truth)
% seperation [um], rel. error., std
%          0    0.0976    0.0813
%     0.4890    0.0634    0.0576
%     0.9780    0.0518    0.0430
%     1.4670    0.0452    0.0390
%     1.9560    0.0441    0.0365
%     2.4450    0.0388    0.0350
%     2.9340    0.0299    0.0320
%     3.4230    0.0222    0.0278
%     4.0750    0.0157    0.0232
%     4.5640    0.0155    0.0197
%     5.0530    0.0159    0.0165

% horzcat(sep_um',ave_sum_abs_val_left_rel_error,std_sum_abs_val_left_rel_error)
% seperation [um], rel. error., std
% 
% ans =
% 
%          0    1.2563    0.1608
%     0.4890    1.3091    0.1569
%     0.9780    1.3481    0.1690
%     1.4670    1.3861    0.1733
%     1.9560    1.4044    0.1795
%     2.4450    1.4129    0.1746
%     2.9340    1.4224    0.1696
%     3.4230    1.4205    0.1653
%     4.0750    1.4271    0.1570
%     4.5640    1.4372    0.1447
%     5.0530    1.4313    0.1478