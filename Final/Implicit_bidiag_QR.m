function B = Implicit_bidiag_QR( B );
    [ m, n ] = size( B );
    
    curm = m;
    while curm > 1
        B(1:curm,1:curm) = Bidiag_Francis_Step( B( 1:curm, 1:curm ) );
        
        
        if abs( B( curm-1,curm ) ) < ... 
            1.0e-14 * ( abs( B( curm-1,curm-1) ) + ...
                        abs( B( curm, curm ) ) )
            curm = curm - 1;
        end
end