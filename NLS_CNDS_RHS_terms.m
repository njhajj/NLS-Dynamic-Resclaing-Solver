%Part of a crank-nicholson solver for the NLS equation in 2 transverse
%dimensions. This function computes various terms that appear on the LHS
%and RHS of the discretized equation.

function[s] = NLS_CNDS_RHS_terms(s)

%linear part of RHS matrix
s.RHS_lin = eye(s.pts) + 1i*s.dxi/2*s.delrho;

%compute the cubic term for predictor or corrector step
    if s.predictor == 1
        s.RHS_nonlin_vec = 1i*s.dxi/2*( 3*abs(s.u).^2.*s.u - 1*abs(s.u_old).^2.*s.u_old );
        s.RHS_a_vec = -1*s.dxi/2*( 3*s.au*s.Drho*(s.rho.*s.u) - 1*s.au_old*s.Drho*(s.rho.*s.u_old) );
    else
        s.RHS_nonlin_vec = 1i*s.dxi/2*(1*abs(s.upred).^2.*s.u + 1*abs(s.u).^2.*s.u_old);
        s.RHS_a_vec = -1*s.dxi/2*( 1*s.aupred*s.Drho*(s.rho.*s.upred) + 1*s.au*s.Drho*(s.rho.*s.u) );
    end

end