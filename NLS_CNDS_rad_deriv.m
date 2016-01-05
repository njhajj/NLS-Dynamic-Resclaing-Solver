%Part of a dynamically scaled crank-nicholson solver for the NLS equation
%in 2 transverse dimensions. This function computes the radial derivates
%that are used. The derivatives are O(dr^4) accurate.

function[s] = NLS_CNDS_rad_deriv(s)

%single radial derivative
Drho = diag(8*ones(s.pts-1,1),1) + diag(-8*ones(s.pts-1,1),-1) + diag(1*ones(s.pts-2,1),-2) + diag(-1*ones(s.pts-2,1),2);

Drho(1,:) = 0; % r=0 BC
Drho(2,2) = 1;

Drho(end-1,end-3) = 0; %r=rmax BC
Drho(end,end-1) = 0;
Drho(end,end-2) = 0;

Drho = Drho/(12*s.drho);
s.Drho = Drho;

%double radial derivative
DDrho = diag(-30*ones(s.pts,1),0) + diag(16*ones(s.pts-1,1),-1) + diag(16*ones(s.pts-1,1),1) + diag(-1*ones(s.pts-2,1),-2) + diag(-1*ones(s.pts-2,1),2); 
DDrho(1,2:end) = 2*DDrho(1,2:end); % r=0 BC
DDrho(2,2) = DDrho(2,2) -1; %BC

Drho(end-1,end-3) = Drho(end-1,end-3)*2; %r=rmax BC
DDrho(end,end-1) = DDrho(end,end-1)*2;
DDrho(end,end-2) = DDrho(end,end-2)*2;

DDrho = DDrho/(12*s.drho^2);

%radial laplacian
s.delrho = zeros(size(DDrho)); %radial laplacian
s.delrho(1,:) = 2*DDrho(1,:); %radial laplacian at r=0 comes from l'hopital's rule
s.delrho(2:end,:) = DDrho(2:end,:) + bsxfun(@times,Drho(2:end,:), 1./s.rho(2:end));

end