clear; clc;
drugs=["Blebbi"];
stiffness=["0.6","1.3","2.6","6","12.7"];
numDrugs=length(drugs);
numStiff=length(stiffness);
%structCell=cell(numStiff,numDrugs);

%for d=[1:numDrugs]
    for s=[1:numStiff]
        if usejava('desktop')
            %[fileSFolders, pathSFolders] = uigetfile('*.mat',['Select ' char(drugs(d)) ' ' char(stiffness(s))]);
            stiffnessFolder=uigetdir('title',['Select ' char(drugs(1)) ' ' char(stiffness(s))]);
        try 
            f=dir(append(stiffnessFolder,'/*.mat'));
            for d=[1:length(f)]
                structCell{s,d}=load(append(stiffnessFolder,'/',f(d).name));
            end
        catch
            disp(['Error :: failed to load file '  stiffnessFolder])
        end
        end
    end

%%SPEED PLOT

t=[0:30]*6;
f=uifigure('Position',[0,0,1000,1000]);
grid=uigridlayout(f,[numDrugs,1]);
first=uipanel('Parent',grid);
first.BorderType='none';
lay=tiledlayout(first,1,numStiff);
for j=[1:numStiff]
    axs(j)=nexttile(lay);
    axs(j).Title.String = ([char(stiffness(j)) 'kPa']);
    axs(j).XLabel.String = ('Time (seconds)');
    axs(j).YLabel.String = ("Speed");
    axs(j).YLim=[0,2000];
    hold(axs(j),'on')
    for i=[1:20]%floor(linspace(1,m,20))
        plot(axs(j),t,structCell{j,1}.speedOut(i,:));
    end
    hold(axs(j),'off')
end
%axs(1).YLim=[0,4000];
%linkaxes(axs,'y')  % scales y to be the same 
ylabel(lay,drugs(1));

for r=[2:numDrugs]
    tempPan=uipanel('Parent',grid);
    tempPan.BorderType='none';
    lay=tiledlayout(tempPan,1,numStiff);
    for j=[1:numStiff]
        axs(j)=nexttile(lay);
        axs(j).XLabel.String = ('Time (seconds)');
        axs(j).YLabel.String = ("Speed");
        axs(j).YLim=[0,2000];
        hold(axs(j),'on')
        for i=[1:20]%floor(linspace(1,m,20))
            plot(axs(j),t,structCell{j,r}.speedOut(i,:));
        end
        hold(axs(j),'off')
    end
%linkaxes(axs,'y')  % scales y to be the same 
    ylabel(lay,drugs(r));


end

%% POWER PLOT
figure()
hold on 
xlabel("Frequency (Hz)");
ylabel("Power")
leg=legend('Location','eastoutside');
title("Power Spectrum of Blebbistatin Speed")
title(leg,'Substrate Stiffness')
bl=structCell(:,1);
for i=[1:length(structCell)]
    allFreq=[];
    for j=[1:length(bl{i}.ps)]
        allFreq(j,:)=rmmissing(bl{i}.ps{j}{2});
    end
    area(bl{i}.ps{1}{1},rescale(mean(allFreq)),'DisplayName',['' char(stiffness(i)) 'kPa'],'FaceAlpha',.5)
end

%% BLEBBI ONLY FIG

t=[0:30]*6;
f2=uifigure('Position',[0,0,1000,200]);
grid2=uigridlayout(f2,[1,1]);
first2=uipanel('Parent',grid2);
first2.BorderType='none';
lay2=tiledlayout(first2,1,numStiff);
for j=[1:numStiff]
    axs2(j)=nexttile(lay2);
    axs2(j).Title.String = ([char(stiffness(j)) 'kPa']);
    axs2(j).XLabel.String = ('Time (seconds)');
    axs2(j).YLabel.String = ("Speed");
    axs2(j).YLim=[0,2000];
    hold(axs2(j),'on')
    
    for k=[1:length(structCell(j,:))]
        if isempty(structCell{j,k})
            continue
        end
        [n,m]=size(structCell{j,k}.speedOut);
        for i=[1:n]%floor(linspace(1,m,20))
            plot(axs2(j),t,structCell{j,k}.speedOut(i,:));
        end
    %hold(axs2(j),'off')
    end
ylabel(lay2,[char(drugs(1)) ])
end
exportapp(f2,'blebbi_speed.pdf')
%% BOX PLOT

%freq vs drug

figure()
boxS=cell(1,numDrugs);
for i=[1:numDrugs]
    stack={};
    for j=[1:length(structCell(:,i))]
        stack=[stack nonzeros(structCell{j,i}.box)];
    end
    boxS{i}=cat(1,stack{:});
end

boxPlotCellArray(boxS,cellstr(drugs))
xlabel("Drug")
ylabel("Frequency (Hz)")

%freq vs stiffness
figure()
boxD=cell(1,numStiff);
for i=[1:numStiff]
%     stack={};
%     for j=[1:length(structCell(i,:))]
%         stack=[stack nonzeros(structCell{i,j}.box)];
%     end
%     boxD{i}=cat(1,stack{:});
    blebbi=structCell{i,1}.box;
    stack=blebbi;
    boxD{i}=cat(1,stack);
end

boxPlotCellArray(boxD,cellstr(stiffness))
xlabel("Stiffness (kPa)")
ylabel("Frequency (Hz)")
