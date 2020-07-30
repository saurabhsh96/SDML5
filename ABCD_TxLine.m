%ABCD Matrix Tx Lines
function A = ABCD_TxLine(Zc, k0, dz)
    A = [cos(k0.*dz) 1j.*Zc.*sin(k0.*dz); 1j*(1./Zc).*sin(k0.*dz) cos(k0.*dz)];
end