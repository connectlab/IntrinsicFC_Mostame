function [ x ] = corr_2D(A, B)
A = A(:); B = B(:);
A( isnan( A ) ) = []; B( isnan( B ) ) = [];
[ acorr, lag ] = xcorr(A - mean( A ), B - mean( B ), 'coeff');
x = acorr( lag == 0 );
end

