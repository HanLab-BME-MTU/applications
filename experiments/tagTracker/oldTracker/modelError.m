function [res, jac]=modelError(parms,match,imgTemp,dataProperties)
%
DEBUG=1;

fSze=size(imgTemp);
for sp=1:length(match)
match(sp).parms=parms(1+(sp-1)*3:3*sp);
end;
%Compute new coords according to model
model=transSpot(match,dataProperties);
%model(s).data=interpn(imgTemp-sl(2).sp(idmap(1).sp{s}).bg,model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3)).*model(s).ratio;
for s =1:length(match)
    model(s).data=interpn(imgTemp-mean(imgTemp(:)),model(s).coords(:,1),model(s).coords(:,2),model(s).coords(:,3),'cubic').*model(s).ratio;
    model(s).residual=-match(s).data+model(s).data;
    %compute gradient
    tImg=zeros(fSze);
    tImg(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)))=model(s).data;
    
    %ZUM TEST:
%     if abs(parms(3))==abs(parms(6))
%         h=findobj(0,'type','figure','name','scan');
%         figure(h);
%         cla;
%         plot(squeeze(imgTemp(model(s).center(1),model(s).center(2),:)));
%         plot(squeeze(tImg(match(s).center(1),match(s).center(2),:)));
%         tImg=zeros(fSze);
%         tImg(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)))=match(s).data;
%         plot(squeeze(tImg(match(s).center(1),match(s).center(2),:)),'r');
%         tImg=zeros(fSze);
%         tImg(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)))=model(s).ratio;
%         plot(squeeze(tImg(match(s).center(1),match(s).center(2),:)),'k-.');
%         pause;
%     end;
    %TEST END
    
    % uncomment following: for check in!!!
    
    [FX,FY,FZ] = gradient(tImg);
    FX=FX(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)));
    FY=FY(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)));
    FZ=FZ(sub2ind(fSze,match(s).coords(:,1),match(s).coords(:,2),match(s).coords(:,3)));
    model(s).iGrad=[FY(:) FX(:) FZ(:)];
end;
%jac=-[FY(:) FX(:) FZ(:)];
res=[];
jac=blkdiag(model.iGrad);
for sp=1:length(match)
    res=[res ; model(sp).residual];
end
h=findobj(0,'Type','figure','Name','res');

% DEBUG
if DEBUG
    global path;
    path=[path;parms];
%     nr1=norm(model(1).residual);
%     nr2=norm(model(2).residual);
%     if ~isempty(h)
%         figure(h);
%         onr=get(h,'UserData');
%         ct=onr(1)+1;
%         plot([onr(1) ct],[onr(2) nr1],'r-');
%         plot([onr(1) ct],[onr(3) nr2],'b-');
%         set(h,'UserData',[ct nr1 nr2]);
         drawnow;
%         if ct>50
%             keyboard;
%         end;
%     else
%         ct=1;
%         h=figure('Name','res');
%         hold on;
%         set(h,'UserData',[ct nr1 nr2]);
%     end;
end;