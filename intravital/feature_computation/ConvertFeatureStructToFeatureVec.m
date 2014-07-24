function [ feature_vec , feature_labels ] = ConvertFeatureStructToFeatureVec( feature_struct , blnOrderFields )

    feature_vec = {};
    feature_labels = {};

    if ~exist( 'blnOrerFields' , 'var' )
    
        blnOrderFields = false;
        
    end
    
    % order structure fields
    if blnOrderFields 
        
        ord_feature_struct = orderfields( feature_struct );
        
    else
        
        ord_feature_struct = feature_struct;
        
    end
    
    % get field name list
    field_name_list = fieldnames( ord_feature_struct );
    
    % convert to feature vec
    for i = 1:numel( field_name_list )
        
        cur_field = getfield( ord_feature_struct , field_name_list{i} );
        
        if isempty( cur_field )
            
             feature_vec{end+1 , 1} = '';
             feature_labels{end+1 , 1} = field_name_list{ i };
            
        elseif isstr( cur_field )

            cur_feature_vec = cur_field;
            cur_feature_labels = field_name_list{ i };
                
        else
            
            if numel( cur_field ) > 1

                cur_field_vals = cur_field(:);
                
                cur_feature_vec = {};
                cur_feature_labels = {};
                
                for v = 1:numel( cur_field_vals )

                    if isstruct( cur_field_vals( v ) )
                        
                        [ vec , labels ] = ConvertFeatureStructToFeatureVec( cur_field_vals( v ) );

                        labels = strcat( field_name_list{ i } , '_' , num2str( v ) , '.' , labels );   
                        
                        cur_feature_vec = [ cur_feature_vec, vec ];
                        cur_feature_labels = [ cur_feature_labels, labels ];            
                            
                     else
                            
                        %cur_feature_vec{ v , 1 } = num2str( cur_field_vals( v ) );
                        cur_feature_vec{ v , 1 } = cur_field_vals( v );
                        
                        cur_feature_labels{ v , 1 } = strcat( field_name_list{ i } , '_' , num2str( v ) );
                        
                     end
                    
                end
                 
            else

                if isstruct( cur_field )
                        
                    [ cur_feature_vec , cur_feature_labels ] = ConvertFeatureStructToFeatureVec( cur_field );

                    cur_feature_labels = strcat( field_name_list{ i } , '.' , cur_feature_labels );

                else                           

                    %cur_feature_vec = num2str( cur_field );
                    cur_feature_vec = cur_field;
                    
                    cur_feature_labels = field_name_list{ i };                

                end

            end      
            
        end
        
        feature_vec = cat(2, feature_vec, (cur_feature_vec(:))' );
        feature_labels = cat(2, feature_labels, (cur_feature_labels(:))');
        
    end    
    
end