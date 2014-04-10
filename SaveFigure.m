function SaveFigure( h , filename , spec_file_format )

    if exist( spec_file_format , 'var' )

        img_file_format = spec_file_format;
        
    else
        
        img_file_format = 'jpeg';
            
    end

    if strcmpi( spec_file_format , 'fig' )
        
        saveas( h, filename, 'fig' ); 
        
    else
        
        style = localGetStyle( h );    
        hgexport( h, filename, style, 'Format' , img_file_format );
        
    end       
    
end

function style = localGetStyle(hfig)

    style = getappdata(hfig,'Exportsetup');
    
    if isempty(style)
      try
        style = hgexport('readstyle','Default');
      catch
        style = hgexport('factorystyle');
      end

    end

end