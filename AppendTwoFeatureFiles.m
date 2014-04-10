function AppendTwoFeatureFiles(file1, file2, outfile)

    numLines = 0;
    
    % copy first file to output
    if exist( 'outfile', 'var' )
        
        fidOut = fopen(outfile, 'w'); 
    
        fidFile1 = fopen(file1, 'r'); 
        strLine = fgetl(fidFile1);
        while ischar(strLine)
           if numLines > 0 
                fprintf(fidOut, '\n%s', strLine); 
           else
                fprintf(fidOut, '%s', strLine);
           end
           numLines = numLines + 1;
           strLine = fgetl(fidFile1);
        end
        fclose(fidFile1);
        
    else
        fidOut = fopen(file1, 'a');
    end
    
    % append second file
    fidFile2 = fopen(file2, 'r'); 
    
    strLine = fgetl(fidFile2); % ignore first line with col headers
    
    strLine = fgetl(fidFile2);
    while ischar(strLine)
       fprintf(fidOut, '\n%s', strLine); 
       numLines = numLines + 1;
       strLine = fgetl(fidFile2);
    end
    fclose(fidFile2);
    fclose(fidOut);
    
end
