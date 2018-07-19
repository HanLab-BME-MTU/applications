function varargout = testTagTracker(test,sourceInfo,targetInfo,movieFrame,constants)
% script to collect the test routines for tag tracker

switch test
    case 1
        % 3D matrix of residuals
        rn=zeros(15,15,8);
        for x=1:15,
            for y=1:15,
                for z=1:8,
                    r=track_lsqnonlinFitFcn(...
                        [-14+(x-1),-14+(y-1),-7+(z-1)]',...
                        sourceInfo,targetInfo,...
                        movieFrame,constants).^2;
                    rn(x,y,z)=mean(r)*100;
                end,
            end,
        end
        imarisShowArray(rn)

    case 2
        mv=zeros(30,30,15,15,3);
        for x=1:15
            [s,t]=extractIntensities(sourceInfo, targetInfo, ...
                movieFrame, constants, [0-(x-1),0,0]');
            f=mv(:,:,:,x,1);
            c=sub2ind([30,30,15],t.coordList(t.goodIdx,1),t.coordList(t.goodIdx,2),t.coordList(t.goodIdx,3));
            f(c)=t.deltaInt;
            mv(:,:,:,x,1)=f;
            tb(x)=t.background;
            mv(:,:,:,x,2)=movieFrame-tb(x);
            f=mv(:,:,:,x,3);
            f(c)=s.intList(s.goodIdx);
            mv(:,:,:,x,3)=f;
        end
        imarisShowArray(mv);
        varargout{1}=mv;
        varargout{2}=tb;

    case 3
        mv=zeros(11,11,9,15);
        for x=1:15
            [s,t]=extractIntensities(sourceInfo, targetInfo, ...
                movieFrame, constants, [0-(x-1),0,0]');
            f=mv(:,:,:,x,1);
            c=sub2ind([11,11,9],...
                t.coordList(t.goodIdx,1)-min(t.coordList(:,1))+1,...
                t.coordList(t.goodIdx,2)-min(t.coordList(:,2))+1,...
                t.coordList(t.goodIdx,3)-min(t.coordList(:,3))+1);
            cc=sub2ind([30,30,15],...
                t.coordList(t.goodIdx,1),...
                t.coordList(t.goodIdx,2),...
                t.coordList(t.goodIdx,3));
            f(c)=t.deltaInt;
            mv(:,:,:,x,1)=f;
            %             tb(x)=t.background;
            %              f=mv(:,:,:,x,2);
            %              f=movieFrame(cc)-tb(x);
            %             mv(:,:,:,x,2)=f;
            %             f=mv(:,:,:,x,3);
            %             f(c)=s.intList(s.goodIdx);
            %             mv(:,:,:,x,3)=f;
        end
        imarisShowArray(mv);
        varargout{1}=mv;

    case 4
        x=[-3:0.1:3]';

        g1=exp(-x.^2/2);
        g11=[g1;zeros(2*length(x),1)];
        [dg,xx]=deal(zeros(length(x),2*length(x)+1));
        for i=1:2*length(x)+1
            xx(:,i)=(i:length(x)+i-1);
            dg(:,i)=g1-g11(xx(:,i));
        end
        cm=jet(2*length(x)+1);
        figure,plot(g11)
        figure,plot(dg);
        colormap(cm);
        figure,set(gca,'NextPlot','add')
        for i=1:2*length(x)+1
            plot(xx(:,i),dg(:,i),'Color',cm(i,:));
        end
        figure,plot(sum(dg.^2,1))

        varargout{1}=dg;

    case 5
        x=[-3:0.1:3]';

        g1=exp(-x.^2/2);
        g1=g1-min(g1);
        g11=[g1;zeros(2*length(x),1)];
        [dg,xx]=deal(zeros(length(x),2*length(x)+1));
        for i=1:2*length(x)+1
            xx(:,i)=(i:length(x)+i-1);
            dg(:,i)=g1-g11(xx(:,i));
        end
        cm=jet(2*length(x)+1);
        figure,plot(g11)
        figure,plot(dg);
        colormap(cm);
        figure,set(gca,'NextPlot','add')
        for i=1:2*length(x)+1
            plot(xx(:,i),dg(:,i),'Color',cm(i,:));
        end
        figure,plot(sum(dg.^2,1))

        varargout{1}=dg;

    case 6
        % 3D matrix of residuals and of convergence
        constants.verbose=-1;
        rn=zeros(51,51,51);
        for x=1:51,
            for y=1:51,
                for z=1:51,
                    xyz = [-1+(x-1)*0.04,-1+(y-1)*0.04,-1+(z-1)*0.04]';
                    r=track_lsqnonlinFitFcn(...
                        xyz,...
                        sourceInfo,targetInfo,...
                        movieFrame,constants).^2;
                    rn(x,y,z)=sum(r);
                end,
            end,
        end
        imarisShowArray(rn)
        varargout{1}=rn;
    case 7
        % coarser representation of residuals
        constants.verbose=-1;
        rn=zeros(161,161,1);
        figure
        for z=1:1,
            for x=1:161,
                for y=1:161,

                    r=track_lsqnonlinFitFcn(...
                        [-8+(x-1)*0.1,-8+(y-1)*0.1,-0+(z-1)*0]',...
                        sourceInfo,targetInfo,...
                        movieFrame,constants).^2;
                    rn(x,y,z)=sum(r);
                    %subplot(2,3,z)

                end,
            end,
            imshow(rn(:,:,z),[])
            colormap(jet);
        end
        [xx,yy]=ndgrid(-8:0.1:8,-8:0.1:8);
        figure,surf(xx,yy,rn)
        %imarisShowArray(rn)

    case 8
        % calculate objective function for different interpolation/readout
        % schemes
%         interpolation = {'*linear',0;...
%             '*linear',1;...
%             '*cubic',0;...
%             '*cubic',1};
        interpolation = {'*cubic',0;...
             '*cubic',1};
        gaussPixOnly = [0,1];
        for int=1:size(interpolation,1)
            constants.interpolation = interpolation(int,:);
            for gauss = 1:2
                constants.gaussPixOnly = gaussPixOnly(gauss);
                testTagTracker(7,sourceInfo,targetInfo,movieFrame,constants);
            end
        end
end
