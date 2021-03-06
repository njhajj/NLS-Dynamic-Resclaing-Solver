%Part of a dynamically scaled crank-nicholson solver for the NLS equation
%in 2 transverse dimensions. This function computes the radial derivates
%that are used. The derivatives are O(dr^2) accurate.

function[s] = NLS_CNDS_rad_deriv_2ndorder(s)

%single radial derivative
Drho = diag(ones(s.pts-1,1),1) + diag(-1*ones(s.pts-1,1),-1);
Drho(1,:) = 0; % r=0 BC
Drho(end,end-3:end) = Drho(end,end-3:end) + [-1,4,-6,4]; % r=rmax boundary includes polynomial extrapolation for external grid point
Drho = Drho/(2*s.drho);
s.Drho = Drho;

%double radial derivative
DDrho = diag(-2*ones(s.pts,1),0) + diag(ones(s.pts-1,1),-1) + diag(ones(s.pts-1,1),1);
DDrho(1,2) = 2*DDrho(1,2); % r=0 BC
DDrho(end,end-3:end) = DDrho(end,end-3:end) + [-1,4,-6,4]; % r=rmax boundary includes polynomial extrapolation for external grid point
DDrho = DDrho/(s.drho^2);
s.DDrho = DDrho;

%radial laplacian
s.delrho = zeros(size(DDrho)); %radial laplacian
s.delrho(1,:) = 2*DDrho(1,:); %radial laplacian at r=0 comes from l'hopital's rule
s.delrho(2:end,:) = DDrho(2:end,:) + bsxfun(@times,Drho(2:end,:), 1./s.rho(2:end));

end