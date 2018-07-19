function save(obj)
save([obj.pathName obj.fileName(1:end-4) '.prv'],'obj');
disp('ROISelector: ROISelector object written to file!');
end