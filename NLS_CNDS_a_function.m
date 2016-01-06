%Part of a crank-nicholson solver for the NLS equation in 2 transverse
%dimensions. This function computes the "a" function used in the solver for
%the current (u), old (u_old), and predicted (upred) value of the nonlinear
%wave u.

%For the psi_mag "a" function I use comes from the paper "Computer
%simulation of wave collapses in the nonlinear Schrodinger equation" where
%it is given in expression 14. With this choice, L(z) = 1/|psi| evaluated
%on axis.

%For the gradpsi_L2 "a" function I use comes from the paper "Focusing
%Singularity of the Cubic Schrodinger Equation" where it is given as
%expression 3.13. For this choice the L function is related to the L2 norm
%of the gradient of the field. In 3.13 I take the arbitrary constant "C" to
%be 1.

function[s] = NLS_CNDS_a_function(s)

    if strcmp(s.Lnorm,'psi_mag') == 1
        
        s.au = -2*imag( (s.DDrho*s.u)./s.u );
        s.au_old = -2*imag( (s.DDrho*s.u_old)./s.u_old );
        s.aupred = -2*imag( (s.DDrho*s.upred)./s.upred );

        s.au = s.au(1);
        s.au_old = s.au_old(1);
        s.aupred = s.aupred(1);
        
    elseif strcmp(s.Lnorm,'gradpsi_L2') == 1
        
        integrand = conj(s.u).^2.*(s.Drho*s.u).^2.*s.rho;
        s.au = imag(trapz(integrand)*s.drho);

        integrand = conj(s.u_old).^2.*(s.Drho*s.u_old).^2.*s.rho;
        s.au_old = imag(trapz(integrand)*s.drho);

        integrand = conj(s.upred).^2.*(s.Drho*s.upred).^2.*s.rho;
        s.aupred = imag(trapz(integrand)*s.drho);
        
    end

end