function popHist(projData)

close all

for i=1:length(projData)

    a=projData(i,1).allTracks; % project number

     segIdx=find(a(:,5)==1);
    fgapIdx=find(a(:,5)==2);
    bgapIdx=find(a(:,5)==3);
    ugapIdx=find(a(:,5)==4);

    figure(i);
    [x1 x2]=hist(a(segIdx,4),100);
    bar(x2,x1,'b')

    hold on

    [x1 x2]=hist(a(fgapIdx,4),100);
    bar(x2,x1,'r')

    [x1 x2]=hist(a(ugapIdx,4),100);
    bar(x2,x1,'g')

    [x1 x2]=hist(a(bgapIdx,4),100);
    bar(x2,x1,'c')

    title('segment(blue), fgap(red), ugap(green), bgap(cyan) velocities')
    xlabel('velocities');
    ylabel('number of segs/gaps');

end


title('segment(red), fgap(blue), bgap(cyan), ugap(green) velocities')


