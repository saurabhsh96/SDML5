%Susceptance Calculations; Infinite
function B = suscpetance_I(omega, dy, dz, m, W)
    eps_0 = 8.854187817e-12;
    constant = omega.*eps_0.*dy./pi;
    term = (abs(sinc((pi.*m.*W./dy)./pi)).^2)./abs(m)...
            .*(1j.*tan(-1j.*pi.*abs(m).*dz./dy));
    B = sum(term).*constant;
end