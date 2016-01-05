%this is for comparing nonlinear phase accumulation of beams with slightly
%different powers. It takes the place of the generic NLS_CNDS_runner
%script.

clear all
close all

% s.pts = 800;
% rho_max = 80;
s.pts = 200;
rho_max = 10;
s.rho = linspace(0,rho_max,s.pts)';
s.drho = s.rho(2) - s.rho(1);
s.dxi = .001/8;

P = 3.3*2; %power in p_crits
u_0 = sqrt(3.77*P);
s.u = u_0*exp(-s.rho.^2);
s.u_old = s.u;
s.upred = s.u;

s.a_xi = 0;
s.L = 1;
s.z = 0;

s = NLS_CNDS_rad_deriv(s);
s = NLS_CNDS_LHS(s);

%create another struct to running with 10% less power

r = s;
P = P/1.1;
u_0 = sqrt(3.77*P);
r.u = u_0*exp(-s.rho.^2);
r.u_old = r.u;
r.upred = r.u;

%%

diff_phase = 0;
m = 1;
n = 1;

while diff_phase < 2*pi
    
    s = NLS_CNDS_stepper(s);
    
    if m == 1
        s.onax_phase(m) = angle(s.psi(1));
    else
        temp = unwrap([s.onax_phase(m-1),angle(s.psi(1))]);
        s.onax_phase(m) = temp(2);
    end
    
    s.onax_psisq(m) = (abs(s.psi(1)).^2);
    
    while r.z(end) < s.z(end)
        r = NLS_CNDS_stepper(r);
        
        if n == 1
            r.onax_phase(n) = angle(r.psi(1));
        else
            temp = unwrap([r.onax_phase(n-1),angle(r.psi(1))]);
            r.onax_phase(n) = temp(2);
        end
        
        r.onax_psisq(n) = (abs(r.psi(1)).^2);
        n = n+1;
    end
    
    %%%%PLOTTING%%%%
    if mod(m,100) == 0
  
        energ = sum(abs(s.u).^2.*s.drho.*s.rho);

        figure(1)
        clf;
        plot(s.r,(abs(s.psi).^2), 'r')
        hold on;
        plot(r.r,(abs(r.psi).^2), 'b')
        title(strcat('intensity vs r at position = ', num2str(m*s.dxi)))

        figure(2)
        clf;
        plot(s.r,log10(abs(s.psi).^2), 'r')
        hold on;
        plot(r.r,log10(abs(r.psi).^2), 'b')
        title(strcat('log of intensity vs r at position = ', num2str(m*s.dxi)))
        
        figure(3)
        clf;
        plot(s.z(2:end),s.onax_phase, 'r')
        hold on;
        plot(r.z(2:end),r.onax_phase, 'b')
        title('on axis phase vs z')
        
        figure(4)
        clf;
        plot(s.r,angle(s.psi), 'r')
        hold on;
        plot(r.r,angle(r.psi), 'b')
        title(strcat('phase vs r, and the energy = ',num2str(energ)))

    end
    
    if m ==100
        ini_energ = energ;
    end
    
    diff_phase = s.onax_phase(end) - r.onax_phase(end);
    m = m+1;
    
end

%%%%finish up propagation of lower power beam so it reaches the same peak
%%%%intensity of the high power beam

while r.onax_psisq(end)<s.onax_psisq(end)
    
    r = NLS_CNDS_stepper(r);
    
    if n == 1
        r.onax_phase(n) = angle(r.psi(1));
    else
        temp = unwrap([r.onax_phase(n-1),angle(r.psi(1))]);
        r.onax_phase(n) = temp(2);
    end
    
    r.onax_psisq(n) = (abs(r.psi(1)).^2);
    n = n+1;
    
        %%%%PLOTTING%%%%
    if mod(n,100) == 0
  
        energ = sum(abs(s.u).^2.*s.drho.*s.rho);

        figure(1)
        clf;
        plot(s.r,(abs(s.psi).^2), 'r')
        hold on;
        plot(r.r,(abs(r.psi).^2), 'b')
        title(strcat('intensity vs r at position = ', num2str(m*s.dxi)))

        figure(2)
        clf;
        plot(s.r,log10(abs(s.psi).^2), 'r')
        hold on;
        plot(r.r,log10(abs(r.psi).^2), 'b')
        title(strcat('log of intensity vs r at position = ', num2str(m*s.dxi)))
        
        figure(3)
        clf;
        plot(s.z(2:end),s.onax_phase, 'r')
        hold on;
        plot(r.z(2:end),r.onax_phase, 'b')
        title('on axis phase vs z')
        
        figure(4)
        clf;
        plot(s.r,angle(s.psi), 'r')
        hold on;
        plot(r.r,angle(r.psi), 'b')
        title(strcat('phase vs r, and the energy = ',num2str(energ)))

    end
    

end