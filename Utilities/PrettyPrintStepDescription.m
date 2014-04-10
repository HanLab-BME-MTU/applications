function PrettyPrintStepDescription( strStepDescription )

    strStar = strStepDescription;
    strStar(:) = '*';
    strStar = [ '****', strStar, '****' ];
    fprintf( '\n\n%s', strStar );
    fprintf( '\n    %s    ', strStepDescription );
    fprintf( '\n%s\n\n', strStar );

end