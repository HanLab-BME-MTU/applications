function circle=convertCircle_xls_to_mat(circle_data_xls)
frameNo=circle_data_xls(:,1);
n=length(frameNo);
for i=1:n
    frame=frameNo(i);
    circle(frame).center=circle_data_xls(i,2:3);
    circle(frame).radius=circle_data_xls(i,4);
end