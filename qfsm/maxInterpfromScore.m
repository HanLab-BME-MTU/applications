function maxV2 = maxInterpfromScore(maxI2,score,vP,vF,mode)
% Sangyoon: I made a change for this refining process to
% use parabola approximation. Once parabola fit is
% too much apart from integer maxV (maxVorg), I
% started to use the fmincon again for more correct refining process.
% parabola approximation
% input:    maxI2       :index for maxV in score
%           score       :score
%           vP,vF       :velocity range
% output:   maxV2       :refined velocity

subv = 1; % radius of subgroup for subscore
maxVorg  = [vP(maxI2(1)) vF(maxI2(2))];

bPolyTracked = 0;
if (maxI2(1)-subv)>=1 && (maxI2(1)+subv)<=size(score,1)...
   && (maxI2(2)-subv)>=1 && (maxI2(2)+subv)<=size(score,2) %|| strcmp(mode,'CDWS')
    subv = 1; % radius of subgroup for subscore
    sub_score = score(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv),...
                        max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));
    % my field of interest
    subvP = vP(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv));
    subvF = vF(max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));

    [subvFG,subvPG]=meshgrid(subvF,subvP);
    subvF1D = reshape(subvFG,[],1);
    subvP1D = reshape(subvPG,[],1);
    sub_score1D = reshape(sub_score,[],1);

    % starting point estimation SH based on discretized maxV (-b/2a =
    % maxVorg(2)) in quadratical expression to avoid the random starting point warning SH
    asp = -0.026; %decided empirically
    bsp = -2*asp*maxVorg(2);
    csp = asp;
    dsp = -2*csp*maxVorg(1);
    esp = -0.5; %arbitrary number
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint', [asp,bsp,csp,dsp,esp]); 
    f = fittype('a*x^2+b*x+c*y^2+d*y+e','independent', {'x', 'y'}, 'dependent', 'z','option',s);
    sf = fit( [subvF1D, subvP1D], sub_score1D, f);

    px = [sf.a sf.b sf.e]; py = [sf.c sf.d sf.e];
    maxV2 = [roots(polyder(py)) roots(polyder(px)) ];
    bPolyTracked = 1;
    
%     % for figure
%     figure
%     hDscore = mesh(subvF,subvP,sub_score);
%     set(hDscore,'EdgeColor',[0 0 0])
%     set(hDscore,'FaceAlpha',0.7)
%     xlabel('displacement in x (pixel)')
%     ylabel('displacement in y (pixel)')
%     zlabel('cross correlation score')
%     zlim([0.987 0.993])
%     set(hDscore,'EdgeColor',[0.6 0.6 0.6])
%     set(hDscore,'LineWidth',2)
%     hold on
%     subfinevF = vF(max(1,maxI2(2)-subv)):0.1:vF(min(size(score,2),maxI2(2)+subv));
%     subfinevP = vP(max(1,maxI2(2)-subv)):0.1:vP(min(size(score,2),maxI2(2)+subv));
%     [subfinevFG,subfinevPG] = meshgrid(subfinevP,subfinevF);
%     hInterp = mesh(subfinevFG,subfinevPG,sf(subfinevFG,subfinevPG));
%     plot3(maxV2(1),maxV2(2),sf(maxV2(1),maxV2(2)),'m.')
%     view(gca,[48.5 10]);
%     set(hInterp,'FaceAlpha',0.3)
    
%     %analytical parabolic fit
%     bap = 0.5*(sub_score(2,3)-sub_score(2,1));
%     cap = 0.5*(sub_score(3,2)-sub_score(1,2));
%     dap = -sub_score(2,2)+0.5*(sub_score(2,3)+sub_score(2,1));
%     eap = -sub_score(2,2)+0.5*(sub_score(3,2)+sub_score(1,2));
%     maxV = [bap/(2*dap),cap/(2*eap)];
%     
%     %analytical cosinusoidal interpolation
%     wx = acos((sub_score(2,2)-sub_score(2,1))/(2*sub_score(2,3)));
%     thetax = atan2((sub_score(2,2)-sub_score(2,1)),(2*sub_score(2,3)*sin(wx)));
%     wy = acos((sub_score(2,2)-sub_score(1,2))/(2*sub_score(3,2)));
%     thetay = atan2((sub_score(2,2)-sub_score(1,2)),(2*sub_score(3,2)*sin(wx)));
%     maxV = [-thetax/wx, -thetay/wy];
end

if (~bPolyTracked || norm(maxVorg-maxV2,2)>1) %&& ~strcmp(mode,'CDWS') %checking for proximity of the fitted point to the original discrete point strcmp(mode, 'accurate') || 
    subv = 4; % expanding region to fit
    sub_score = score(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv),...
                        max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));
    % my field of interest
    subvP = vP(max(1,maxI2(1)-subv):min(size(score,1),maxI2(1)+subv));
    subvF = vF(max(1,maxI2(2)-subv):min(size(score,2),maxI2(2)+subv));

    sp   = csape({subvP,subvF},sub_score);
    dsp1 = fnder(sp,[1,0]);
    dsp2 = fnder(sp,[0,1]);
    options = optimset('Algorithm','interior-point'); % this generates an warning
%             options = optimset('Algorithm','sqp');% this too

    maxV2 = fmincon(@vFun,maxVorg,[],[],[],[], ...
        [max(subvP(1),maxVorg(1)-2) max(subvF(1),maxVorg(2)-2)], ...
        [min(subvP(end),maxVorg(1)+2), min(subvF(end),maxVorg(2)+2)],[], ...
        options,sp,dsp1,dsp2);            
end
