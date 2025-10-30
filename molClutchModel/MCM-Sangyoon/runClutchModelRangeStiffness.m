function [mf_G,v_G,mf_GErr,v_GErr] = ...
    runClutchModelRangeStiffness(nExp,ksub,nm,fm1,vu,nc,dint1,dint2,kont1,...
           kont2,kof1,kof2,kc,konv,pt,mr,intadd,ion,v_actin,dActin,tTotal)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

mfGroup = cell(1,nExp);
mvGroup = cell(1,nExp);
% mnb1GroupCK1 = cell(1,nExp);
% mnb2GroupCK1 = cell(1,nExp);
% mdint1GroupCK1 = cell(1,nExp);
% mdint2GroupCK1 = cell(1,nExp);
numKsub = length(ksub);
mf = zeros(numKsub,1);
mv = zeros(numKsub,1);

for p=1:nExp
    for ii=1:numKsub
%        [mfi,mvi,mnb1i,mnb2i,mdint1i,mdint2i] = ...
       [mfi,mvi] = clutchModelNascentAdhesion(nm,fm1,vu,nc,dint1,dint2,kont1,...
           kont2,kof1,kof2,kc,ksub(ii),konv,pt,mr,intadd,ion,v_actin,dActin,tTotal);
        mf(ii) = mfi;
        mv(ii) = mvi;
%         mnb1(ii) = mnb1i;
%         mnb2(ii) = mnb2i;
%         mdint1(ii) = mdint1i;
%         mdint2(ii) = mdint2i;
        disp([num2str(100*ii/numKsub) '% done...'])
    end
    mfGroup{p} = mf;
    mvGroup{p} = mv;
%     mnb1GroupCK1{p} = mnb1;
%     mnb2GroupCK1{p} = mnb2;
%     mdint1GroupCK1{p} = mdint1;
%     mdint2GroupCK1{p} = mdint2;
end
v_G = mean(cell2mat(mvGroup),2); %mv;
mf_G = mean(cell2mat(mfGroup),2); %mf;
% mdint1_ck666 = mean(cell2mat(mdint1GroupCK1),2); %mdint1;
mf_GErr = std(cell2mat(mfGroup),0,2);
v_GErr = std(cell2mat(mvGroup),0,2);

end