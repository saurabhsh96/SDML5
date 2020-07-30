%Susceptance Calculations
function B = suscpetance(omega, dy, dx, m, W)
    eps_0 = 8.854187817e-12;
    constant = omega.*eps_0.*dy./pi;
    term = 0;
    for ind = 1:length(m)
        if (m(ind) == 0)
            %disp('here');
            continue;
        end
        term = term + (abs(sinc((pi*m(ind).*W./dy)./pi)).^2)./abs(m(ind));
    end
    B = term.*constant;
end