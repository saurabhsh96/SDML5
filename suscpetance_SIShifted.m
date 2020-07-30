%Susceptance Calculations; Semi-Infinite
function B = suscpetance_SIShifted(omega, dy, dz, m, W, sy)
%     eps_0 = 8.854187817e-12;
    c = 3e8;
    zeta0 = 120.*pi;
    constant = omega./c.*dy./pi./zeta0;
    term = (abs(sinc((pi*m.*W./dy)./pi)).^2)./abs(m)...
        .*(1./2-1j./2.*cot(-2j.*pi.*abs(m).*dz./dy)...
            + 1j./2.*exp(2j.*pi.*m.*sy./dy).*csc(-2j.*pi.*abs(m).*dz./dy));
    B = sum(term).*constant;
end