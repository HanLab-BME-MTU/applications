function connectivity=getConnectivity(point1, point2, labeledImage)

%point1 and point2 have to be x and y coordinate matrices of the form
%[Xpoint Ypoint; Xnei1 Ynei1; Xnei2 Ynei2,Xnei3 Ynei3]

% Author: Santiago Costantino 
% santiago.costantino@umontreal.ca


connectivity=zeros(1,size(point1, 1)+size(point2, 1)-2);
counter=1;
for itPoint1=1:size(point1, 1)-1
    for itPoint2=1:size(point2, 1)-1
        counter=counter+1;
        if labeledImage(point1(itPoint1+1,2),point1(itPoint1+1,1))==...
            labeledImage(point2(itPoint2+1,2),point2(itPoint2+1,1))
            connectivity(counter)=labeledImage(point2(itPoint2+1,2),point2(itPoint2+1,1));
        else
            % Got to take into account the cases in which two joints are
            % too close together with no segment in between. In that case
            % the points have only three elements. The point and two
            % neighbour pixels
            if and(size(point1, 1)==3, size(point2, 1)==3)
                distance=sqrt((point1(1,1)-point2(1,1))^2+(point1(1,2)-point2(1,2))^2);
                if distance<7
                    connectivity(counter)=1;
                end
            else
                connectivity(counter)=0;
            end
        end
    end
end

% for itObjects=1:labelNumbers
%     testImage=zeros(size(labeledImage));
%     testImage(labeledImage==itObjects)=1;
%     testImage(point1(1, 2),point1(1, 1))=1;
%     testImage(point2(1, 2),point2(1, 1))=1;
%     [testImageLabeled, thisNumber]=bwlabel(testImage, 8);
%     if thisNumber==1
%         connectivity=1;
%         return
%     else
%         connectivity=0;
%     end
% end
    
    