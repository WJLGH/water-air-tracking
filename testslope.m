Hkx = H0.*1i.*Kx;
Hky = H0.*1i.*Ky;
Zx = real( ifft2(Hkx) .* Grid_Sign ) ;
Zy = real( ifft2(Hky) .* Grid_Sign ) ;
n = sqrt(Zx.^2+Zy.^2+1);
Zx = Zx./n;
Zy = Zy./n;