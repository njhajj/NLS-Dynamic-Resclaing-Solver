%Part of a crank-nicholson solver for the NLS equation in 2 transverse
%dimensions. This function steps forward the solution

function[s] = NLS_CNDS_stepper(s)

% pts = length(s.r);

s.predictor = 1; %switch to predictor step
s = NLS_CNDS_a_function(s); %compute the "a" function
s = NLS_CNDS_RHS_terms(s); %compute terms used in the inversion

%predictor step
RHS_vec = s.RHS_lin*s.u + s.RHS_nonlin_vec + s.RHS_a_vec;
s.upred = s.inv_LHS*RHS_vec; %1/4 of the computation time

s.predictor = 0; %switch to corrector step
s = NLS_CNDS_a_function(s); %compute the "a" function
s = NLS_CNDS_RHS_terms(s); %compute terms used in the inversion

%corrector step
RHS_vec = s.RHS_lin*s.u + s.RHS_nonlin_vec + s.RHS_a_vec;
u_out = s.inv_LHS*RHS_vec;

%update fields
s.u_old = s.u;
s.u = u_out;
% s.u(end) = 0;
s.a_xi(end+1) = (1/2)*(s.au +s.aupred);
s = NLS_CNDS_inverter(s);


end