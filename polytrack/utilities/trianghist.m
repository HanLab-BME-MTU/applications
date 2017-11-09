function [distances,avaragedist] = trianghist(coor,counter)
%Make a Delauneytriangulation of a set of coordinates and
%calculate the distances between points without redundancies

tri=delaunay(coor(:,1),coor(:,2));
%returns the indexes of the coordinates, that make up a triangle together

TRI=tri;
TRI(:,4)=tri(:,1);

%match pairs of points, describing one line (no longer three points, describing a
%triangle)and chuck out the lines(between two points) that occure twice 
uni=cat(1,unique(TRI(:,1:2),'rows'),unique(TRI(:,2:3),'rows'),unique(TRI(:,3:4),'rows'));
distances=zeros(length(uni),1);

%calculate the distances
for i=1:length(uni)   
        distances(i)=sqrt((coor(uni(i,1),1)-coor(uni(i,2),1))^2+(coor(uni(i,1),2)-coor(uni(i,2),2))^2);  
end

avaragedist=sum(distances)/length(distances);

figure, hist(distances,[1:10:900]), title(num2str(counter));
axis([0 900 0 35]);
hold on
text(700,22,num2str(avaragedist));
text(700,18,num2str(length(distances)));
hold off