function trihistgroups(infos)
%this programm takes the infos from master and calculates the 
%distances between the cells belonging to the same group(this
%information is within infos)


%for every picture there is one field of the type info
for h=1:size(infos,2)
    distances=[0];
    clear members;
    clear comemb ;
  % the third row of infos.info says, to which group this nuc belongs
    groups=unique(infos(h).info(:,3));
    
    for k=1:length(groups)
         %find the members and the coord of these members for every group
         members=find(infos(h).info(:,3)==groups(k));
         comemb=infos(h).info(members,1);
         comemb(:,2)=infos(h).info(members,2);
         
         if size(comemb,1)>2
             %more than two cells in the group; triang possible
           [distancestemp,gulp]=trianghist(comemb,k);
           %add the distances between the cells of this group to the
           %distances of the cells in other groups(same image)
           distances=cat(1,distances,distancestemp);
       end
     end
     
%make a histogramm of the distances         
figure, hist(distances,[1:10:900]), title(num2str(h));
axis([0 600 0 length(distances)/5]);
hold on

%how many lines between a pair of cells are there for this picture
text(500,20,'total count of distances between a pair of cells');
text(500,18,num2str(length(distances)));
hold off

%calculate some stuff for this picture and store the data
y(h,1) = skewness(distances);
avaragedist(h,1)=sum(distances)/length(distances);    
kurt(h,1)=kurtosis(distances);
end

%plot the accumulated calculated stuff for all pictures
figure, title('Only nucloi within groups: skewness/blue   avaragedist/green   kurtosis/red')
hold on
plot(y)
plot(avaragedist/100,'g')
plot(kurt/4,'r')
hold off

