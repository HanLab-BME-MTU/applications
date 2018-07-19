function fixEpsFile(filename)

[p,b,n,ext] = getFilenameBody(filename);

if ~strcmp(ext,'.eps')
    error([filename ' is not a EPS file.']);
end

tmpname = [p filesep b n '_TMP' ext];

fid_in = fopen(filename,'r');
fid_out = fopen(tmpname,'w');

tline = fgetl(fid_in);
while ischar(tline)
    stridx = strfind(tline,'\302\265');
    if ~isempty(stridx)
        tline = [tline(1:stridx) '265' tline(stridx+8:end)];
    end
    fprintf(fid_out, '%s\n',tline);
    tline = fgetl(fid_in);
end

fclose(fid_in);
fclose(fid_out);
movefile(tmpname,filename);
