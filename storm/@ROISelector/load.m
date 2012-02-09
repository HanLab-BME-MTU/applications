function roiSel = load(fullPath)
tmp = load([fullPath(1:end-4) '.prv'],'-mat');
roiSel = tmp.obj;
disp('ROISelector: ROISelector object read from file!');
end