%Function to convert ABCD to S
function S = ABCDtoS(Fmat, Z0)
    denom = Fmat(1,1) + Fmat(1,2)./Z0 + Fmat(2,1).*Z0 + Fmat(2,2);   
    S(1, 1) = (Fmat(1,1) + Fmat(1,2)./Z0 - Fmat(2,1).*Z0 - Fmat(2,2))./denom;
    S(1, 2) = 2.*(Fmat(1,1).*Fmat(2,2) - Fmat(1,2).*Fmat(2,1))./denom;
    S(2, 1) = 2./denom;
    S(2, 2) = (-Fmat(1,1) + Fmat(1,2)./Z0 - Fmat(2,1).*Z0 + Fmat(2,2))./denom;
end