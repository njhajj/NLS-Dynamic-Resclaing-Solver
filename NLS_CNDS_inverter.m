%Part of a crank-nicholson solver for the NLS equation in 2 transverse
%dimensions. This function transforms the dynamically rescaled solution
%u(xi,rho) to the solution of the nondimensional NLS in its original
%coordinates psi(r,z).

function[s] = NLS_CNDS_inverter(s)
    
    %upgrade L
    xi_reg = [s.a_xi(end-1),s.a_xi(end)];
    s.L(end+1) = s.L(end)*exp(-trapz(xi_reg)*s.dxi);
    
    %upgrade variables for the "original" non-dimensional 2D NLS problem
    s.r = s.rho*s.L(end);    
    s.psi = (1/s.L(end))*s.u;
    s.z(end+1) = s.z(end) + s.L(end)^2*s.dxi;

end