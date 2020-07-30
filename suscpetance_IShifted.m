%Susceptance Calculations; Infinite
function B = suscpetance_IShifted(omega, dy, dz, m, W, sy)
    %eps_0 = 8.854187817e-12;
    zeta0 = 120.*pi;
    c = 3e8;
    constant = omega./c.*dy./pi./zeta0;
    term = (abs(sinc((pi.*m.*W./dy)./pi)).^2)./abs(m)...
            .*(-1j.*cot(-2j.*pi.*abs(m).*dz./dy)...
            + 1j.*exp(2j.*pi.*m.*sy./dy).*csc(-2j.*pi.*abs(m).*dz./dy));
    B = sum(term).*constant;
end