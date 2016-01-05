%Part of a crank-nicholson solver for the NLS equation in 2 transverse
%dimensions. This function computes the "a" function used in the solver for
%the current (u), old (u_old), and predicted (upred) value of the nonlinear
%wave u.

%The particular choice of "a" function I use comes from the paper "Focusing
%Singularity of the Cubic Schrodinger Equation" where it is given as
%expression 3.13.

function[s] = NLS_CNDS_a_function(s)

integrand = conj(s.u).^2.*(s.Drho*s.u).^2.*s.rho;
s.au = imag(trapz(integrand)*s.drho);
% s.au = imag(sum(integrand)*s.drho);

integrand = conj(s.u_old).^2.*(s.Drho*s.u_old).^2.*s.rho;
s.au_old = imag(trapz(integrand)*s.drho);
% s.au_old = imag(sum(integrand)*s.drho);

integrand = conj(s.upred).^2.*(s.Drho*s.upred).^2.*s.rho;
s.aupred = imag(trapz(integrand)*s.drho);
% s.aupred = imag(sum(integrand)*s.drho);

end