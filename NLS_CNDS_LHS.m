%Part of a crank-nicholson solver for the NLS equation in 2 transverse
%dimensions. This function computes various terms that appear on the LHS
%and RHS of the discretized equation.

function[s] = NLS_CNDS_LHS(s)

    %LHS matrix
    s.LHS = sparse(eye(s.pts) - 1i*s.dxi/2*s.delrho);
    
%     %r=rmax boundary conditions
%     s.LHS(end,:) = 0;
%     s.LHS(end-1,:) = 0;
%     s.LHS(end,end) = 1;
%     s.LHS(end-1,end-1) = 1;
    
    %invert the LHS
    s.inv_LHS = full(inv(s.LHS));

end