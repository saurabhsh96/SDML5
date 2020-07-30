%Lecture 5, Question 1 For single layer
clear;
close all;
%% Defining Inputs

freq = 2e9:0.5e9:10e9;
c = 3e8;
lam = c./freq;
k0 = 2*pi./lam;

fr0 = 10e9;
lam0 = c./fr0;
W = 0.01.*lam0;
dx = 0.2.*lam0;
dy = dx;

drad = pi/180;
th = [0, 60].*drad;
phi = 0;

eps_0 = 8.854187817e-12;
mu_0 = 1.2566370614e-6;

zeta0 = 120*pi;
%% Susceptance

%m = [-4:-1, 1:5];
m = -10:1:10;
B = zeros(1, length(freq));

for ind = 1:length(freq)
    omega = 2*pi*freq(ind);
    B(ind) = suscpetance(omega, dy, dx, m, W);
end

%% S11 Calculations

%For different incident angles
for ind = 1:length(th)
    ZTM = -1j./B;
    ZTE = ZTM.*(1./(1-((sin(th(ind)).^2)./2)));
    
    kx = k0.*sin(th(ind)).*cos(phi);
    ky = k0.*sin(th(ind)).*sin(phi);
    kRho = sqrt(kx.^2 + ky.^2);
    kz0 = -1j*sqrt(-((k0.^2)-(kRho.^2)));
    
    Z0TE = (zeta0.*k0)./kz0;
    Z0TM = (zeta0.*kz0)./k0;
    
    ZTMReq = (ZTM.*Z0TM)./(ZTM + Z0TM);
    ZTEReq = (ZTE.*Z0TE)./(ZTE + Z0TE);
    
    S11TM = (Z0TM - ZTMReq)./(Z0TM + ZTMReq);
    S11TE = (Z0TE - ZTEReq)./(Z0TE + ZTEReq);
    
    S12TM = 1 - abs(S11TM).^2;
    S12TE = 1 - abs(S11TE).^2;
    
    %Plotting
    figure();
    plot(freq./10^9, abs(S11TM).^2, 'LineWidth', 1.5, 'DisplayName', 'S11(TM)'); hold on;
    plot(freq./10^9, abs(S11TE).^2, 'LineWidth', 1.5, 'DisplayName', 'S11(TE)');
    plot(freq./10^9, abs(S12TM), 'LineWidth', 1.5, 'DisplayName', 'S12(TE)');
    plot(freq./10^9, abs(S12TE), 'LineWidth', 1.5, 'DisplayName', 'S12(TE)');
    title(['|S11|^2 and |S12|^2 vs. Frequency at, theta =', num2str(th(ind))]);
    xlabel('Frequency (in GHz)');
    ylabel('|S11|^2 and |S12|^2');
    legend show;
    grid on;
    hold off;
end


