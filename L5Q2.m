%Lecture 5, Question 1 For single layer
clear;
close all;
%% Defining Inputs

%EM Specs
freq = 2e9:0.5e9:10e9;
c = 3e8;

fr0 = 10e9;
lam0 = c./fr0;
W = 0.01.*lam0;
dx = 0.2.*lam0;
dy = dx;
dz = W;

%Number of layers
N = 5;

%Space
drad = pi/180;
th = [0, 60].*drad;
phi = 0;

%Impedance
eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;

zeta0 = 120*pi;

m = [-10:-1, 1:10];

%% Calculating ABCD Matrices

for ind = 1:length(th)
    for indF = 1:length(freq)
        %Susceptance
        omega = 2*pi*freq(indF);
        B_SI = suscpetance_SI(omega, dy, dz, m, W);
        B_I = suscpetance_I(omega, dy, dz, m, W);

        lam = c./freq(indF);
        k0 = 2*pi./lam;
        
        %Propagation constant
        kx = k0.*sin(th(ind)).*cos(phi);
        ky = k0.*sin(th(ind)).*sin(phi);
        kRho = sqrt(kx.^2 + ky.^2);
        kz0 = -1j*sqrt(-((k0.^2)-(kRho.^2)));

        %Tx Line Impedance
        Z0TE = (zeta0.*k0)./kz0;
        Z0TM = (zeta0.*kz0)./k0;

        %Edge impedance
        ZTMedge = -1j./B_SI;
        ZTEedge = ZTMedge.*(1./(1-((sin(th(ind)).^2)./2)));

        %Middle impedance
        ZTMinf = -1j./B_I;
        ZTEinf = ZTMinf.*(1./(1-((sin(th(ind)).^2)./2)));

        %Matrices
        matEdgeTM = ABCD_Z(ZTMedge); 
        matEdgeTE = ABCD_Z(ZTEedge);

        matBetTM = ABCD_Z(ZTMinf);
        matBetTE = ABCD_Z(ZTEinf);

        matTxTM = ABCD_TxLine(Z0TM, k0, dz);
        matTxTE = ABCD_TxLine(Z0TE, k0, dz);

        %Final Matrix
        FmatTM = (matEdgeTM^2)*(matBetTM^(N-2))*(matTxTM^(N-1));
        FmatTE = (matEdgeTE^2)*(matBetTE^(N-2))*(matTxTE^(N-1));

        Stm = ABCDtoS(FmatTM, Z0TM);
        Ste = ABCDtoS(FmatTE, Z0TE);
        
        S11tm(indF) = Stm(1,1);
        S11te(indF) = Ste(1,1);
        S12tm(indF) = Stm(1,2);
        S12te(indF) = Ste(1,2);
    end
    %Plotting
    figure();
    plot(freq./10^9, abs(S11tm).^2, 'LineWidth', 1.5, 'DisplayName', 'S11(TM)'); hold on;
    plot(freq./10^9, abs(S11te).^2, 'LineWidth', 1.5, 'DisplayName', 'S11(TE)');
    plot(freq./10^9, abs(S12tm).^2, 'LineWidth', 1.5, 'DisplayName', 'S12(TE)');
    plot(freq./10^9, abs(S12te).^2, 'LineWidth', 1.5, 'DisplayName', 'S12(TE)');
    title(['|S11|^2 and |S12|^2 vs. Frequency at, theta =', num2str(th(ind))]);
    xlabel('Frequency (in GHz)');
    ylabel('|S11|^2 and |S12|^2');
    legend show;
    grid on;
    hold off;
end