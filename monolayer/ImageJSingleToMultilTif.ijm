dir = getDirectory("Choose a Directory "); 
saveDir = getDirectory("Choose a Directory "); 
listFiles(dir);
  
  function listFiles(dir) {
     list = getFileList(dir);
     for (i=0; i<list.length; i++) {
        if (endsWith(list[i], "/")){
           listFiles(""+dir+list[i]);}
        else {
for (i = 0; i<list.length; i++){
	if (endsWith(list[i], "p000001t00000001z001c04.tif")){
run("Image Sequence...", "open="+dir+list[i]+" sort");
title = getInfo("slice.label");
string2 = replace(title, "_p000001t00000001z001c0", "");
string3 = replace(string2,"-","_");
saveAs("tif", saveDir+File.separator+string3);
close(string3+".tif");}
}}}}
