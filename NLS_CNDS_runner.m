clear all
close all

s.pts = 400;
rho_max = 40;
s.rho = linspace(0,rho_max,s.pts)';
s.drho = s.rho(2) - s.rho(1);
s.dxi = .001/8;

P = 1*2; %power in p_crits
u_0 = sqrt(3.77*P);
u_0 = 4;
s.u = u_0*exp(-s.rho.^2);
% s.u = 4./(1.+s.rho).^2;
% s.u(end) = 0;
s.u_old = s.u;
s.upred = s.u;

s.a_xi = 0;
s.L = 1;
s.z = 0;

s = NLS_CNDS_rad_deriv(s);
s = NLS_CNDS_LHS(s);

%%
onax_psisq = 0;
m = 1;
while onax_psisq(end)<60000
% for m = 1:10*round(1/s.dxi)
    s = NLS_CNDS_stepper(s);
    onax_phase(m) = angle(s.psi(1));
    onax_psisq(m) = (abs(s.psi(1)).^2);
    if mod(m,100) == 0
        
        energ = sum(abs(s.u).^2.*s.drho.*s.rho);
        figure(1)
        clf;
        plot(s.r,(abs(s.psi).^2), 'r')
        title(strcat('intensity vs r at position = ', num2str(m*s.dxi)))

        figure(2)
        clf;
        plot(s.r,log10(abs(s.psi).^2), 'r')
        title(strcat('log of intensity vs r at position = ', num2str(m*s.dxi)))
        
        figure(3)
        clf;
        plot(s.z(2:end),onax_phase)
        title('on axis phase vs z')
        
        figure(4)
        clf;
        plot(s.r,angle(s.psi))
        title(strcat('phase vs r, and the energy = ',num2str(energ)))
        
%         pause;
        
        
    end
    m = m+1;

end
