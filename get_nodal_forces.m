function [ Pr ] = get_nodal_forces(elems, u, L, T, fy, ey, n, ndof, rests)
% This function calculates the resisting nodal forces based on the stress
% in the elements.
%Author: Nicole Widmer
%Last update: 23/10/2021

    Pr = zeros(size(u));

    for m = 1:n     %for every bar

        % Node 1 of element
        N1=elems(m,1);

        % Node 2 of element
        N2=elems(m,2);

        % find indices of DOFs
        i1=(N1*2-1); i2=N1*2; % indeces of DOFs of Node 1
        i3=(N2*2-1); i4=N2*2; % indeces of DOFs of Node 2
        ii=[i1 i2 i3 i4];

        % Element properties
        A = elems(m, 3);        %area [m2]
        E = elems(m, 4);        %E-modulus [Pa]
        
        % local displacement of bar
        u_loc = T(:,:,m)*u(ii);

        % strain in bar
        eps = [-1 0 1 0] / L(m) * u_loc;

        %calculate stress in bar based on strain and yield limit
        if abs(eps) >= ey
            p = fy * A + E*A*(abs(eps)-ey);
            p = sign(eps) * p;
        else
            p = E*A*eps;
        end 

        %calculate global nodal force
        P = T(:,:,m)' * [-1 0 1 0]' * p;
        Pr(ii) = Pr(ii) + P;
    end

    % set nodal forces of restrained DOF to zero
    for j = 1:ndof
        if rests(j) == 1
            Pr(j) = 0;
        end
    end

end

