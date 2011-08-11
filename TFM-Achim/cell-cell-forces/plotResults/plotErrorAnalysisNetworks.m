function []=plotErrorAnalysisNetworks(groupedClusters)
cutOffMag=0;
cutOffRelErr=Inf;

%goodEdgeSet=findEdges(groupedClusters,'kPa',[8 35],'relErrF',Inf,'errs',0);
goodEdgeSet=findEdges(groupedClusters,'relErrF',Inf,'errs',0);
[f1_vals,f2_vals,fc1_vals]=collectEdgeValues(groupedClusters,goodEdgeSet,'f1','f2','fc1','fc2');
f1_mag = sqrt(sum(f1_vals.^2,2));
f2_mag = sqrt(sum(f2_vals.^2,2));
fc_mag = sqrt(sum(fc1_vals.^2,2));
fn_vals= 0.5*(f1_vals-f2_vals);
fn_mag = sqrt(sum(fn_vals.^2,2));

% remove the edges that can not be covered by the network analysis:
checkVec=isnan(fn_mag);
f1_mag(checkVec)=[];
f2_mag(checkVec)=[];
fc_mag(checkVec)=[];
fn_mag(checkVec)=[];
fn_vals(checkVec,:)=[];
f1_vals(checkVec,:)=[];
f2_vals(checkVec,:)=[];
fc1_vals(checkVec,:)=[];

figure()
glbMin=min([f1_mag;f2_mag]);
glbMax=max([f1_mag;f2_mag]);
% plot(f1_mag,f2_mag,'.b');
plot(fn_mag,f1_mag,'.b')
hold on
plot(fn_mag,f2_mag,'.b')
% these are the abs erors:
% plot(fn_mag,sqrt(sum(( f1_vals-fn_vals).^2,2)),'.b')
% hold on;
% plot(fn_mag,sqrt(sum((-f2_vals-fn_vals).^2,2)),'.b')
plot(fc_mag,f2_mag,'.r');
plot(fc_mag,f1_mag,'.r');
plot([0 glbMax],[0 glbMax],'--k')
title('magnitude of fi vs <f_{net}> [b] or fc [r]')
xlim([glbMin glbMax])
ylim([glbMin glbMax])
xlabel('<f_{net}> or fc [nN]')
ylabel('f1 or f2 [nN]')
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off


checkVec=fn_mag<cutOffMag;
f1_mag(checkVec)=[];
f2_mag(checkVec)=[];
fn_mag(checkVec)=[];
fc_mag(checkVec)=[];
f1_vals(checkVec,:)=[];
f2_vals(checkVec,:)=[];
fc1_vals(checkVec,:)=[];
fn_vals(checkVec,:)=[];

alpha_f12_to_fn= vertcat(acosd(dot(f1_vals,fn_vals,2)./(f1_mag.*fn_mag)),acosd(dot(-f2_vals,fn_vals,2)./(f2_mag.*fn_mag)));
% alpha_f1_f2= acos(dot(f1_vals,f2_vals,2)./(f1_mag*f2_mag));
alpha_f12_to_fc= vertcat(acosd(dot(f1_vals,fc1_vals,2)./(f1_mag.*fc_mag)),acosd(dot(-f2_vals,fc1_vals,2)./(f2_mag.*fc_mag)));

display(['fn: ',num2str(sum(alpha_f12_to_fn>90)),'(= ',num2str(sum(alpha_f12_to_fn>90)/numel(isnan(alpha_f12_to_fn))),') points are >90deg']);
display(['fc: ',num2str(sum(alpha_f12_to_fc>90)),'(= ',num2str(sum(alpha_f12_to_fc>90)/numel(isnan(alpha_f12_to_fc))),') points are >90deg']);


figure()
hist(alpha_f12_to_fn,250,'b')
title('Angular histogram fn')
xlabel('angle')
ylabel('counts')
xlim([0 180])
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off

figure()
hist(alpha_f12_to_fc,250,'r')
title('Angular histogram fc')
xlabel('angle')
ylabel('counts')
xlim([0 180])
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off

%close all;

% append f1,f2 and doublicate the other measurements:

fi_vals_all = vertcat(f1_vals,-f2_vals);
fi_mag_all  = vertcat(f1_mag , f2_mag);

fn_vals_all = vertcat(fn_vals, fn_vals);
fn_mag_all  = vertcat(fn_mag , fn_mag);

fc_vals_all= vertcat(fc1_vals,fc1_vals);
fc_mag_all  = vertcat(fc_mag  ,fc_mag);

% Gaudenz suggestion:
% rel_err_fn = 2*(f1_mag-f2_mag)./(f1_mag+f2_mag);
% better suggestion:
% rel_err_fn = (fi_mag_all-fn_mag_all)./fn_mag_all;
% my suggestion:
rel_err_fn = sqrt(sum((fi_vals_all-fn_vals_all).^2,2))./fn_mag_all;
figure()
hist(rel_err_fn,2000);
%xlim([0 1])
title('rel. error: (fi_{mag}-<f_{net}>)./(<f_{net}>)')
xlabel('rel. error')
ylabel('counts')
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off

rel_err_fn(abs(rel_err_fn)>cutOffRelErr)=[];
display('rel. error fn:')
display(['mean   +-      std: ',num2str(  mean(rel_err_fn),'%0.4f'),' +- ',num2str(       std(rel_err_fn),'%0.4f')])
display(['median +- 1.48*MAD: ',num2str(median(rel_err_fn),'%0.4f'),' +- ',num2str(1.4826*mad(rel_err_fn),'%0.4f')])



% Old analysis:
% rel_err_fc_fn = (fn_mag_all-fc_mag_all)./(fn_mag_all+fc_mag_all)/2;
% better suggestion:
% rel_err_fc_fn = (fi_mag_all-fc_mag_all)./fc_mag_all;
% my suggestion:
rel_err_fc_fn = sqrt(sum((fi_vals_all-fc_vals_all).^2,2))./fc_mag_all;

figure()
% checked with part7
% hist(rel_err_fn ,500);
hist(rel_err_fc_fn ,2000);
title('rel. error: 2*(fn_{mag}-fc_{mag})./(fn_{mag}+fc_{mag})')
%title('rel. error: (fi_{mag}-fc_{mag})./(fc_{mag})')
xlabel('rel. error')
ylabel('counts')
%xlim([0 1])
ylim([0 550])
box on
set(gca,'LineWidth',2,'FontSize',20)
hold off

rel_err_fc_fn(abs(rel_err_fc_fn)>cutOffRelErr)=[];
display('rel. error fc:')
display(['mean   +-      std: ',num2str(  mean(rel_err_fc_fn),'%0.4f'),' +- ',num2str(       std(rel_err_fc_fn),'%0.4f')])
display(['median +- 1.48*MAD: ',num2str(median(rel_err_fc_fn),'%0.4f'),' +- ',num2str(1.4826*mad(rel_err_fc_fn),'%0.4f')])
