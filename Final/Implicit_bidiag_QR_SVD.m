function [ U, B, V ] = Implicit_bidiag_QR_SVD( U, B, V );
    [ m, n ] = size( B );
    
    curm = m;
    while curm > 1
        [ U, B(1:curm,1:curm), V ] = ...
            Bidiag_Francis_Step_Update_U_V( U, B( 1:curm, 1:curm ), V );
        
        if abs( B( curm-1,curm ) ) < ... 
            1.0e-14 * ( abs( B( curm-1,curm-1) ) + ...
                        abs( B( curm, curm ) ) )
            curm = curm - 1;
        end
end