function [epsilon] = get_strain(elems, u, L, T, element)
% This function calculates the engineering strain in a bar element based on the displacements of the nodes
% and the initial length of the element.
%Author: Nicole Widmer
%Last update: 23/10/2021

    node_element=[elems(element,1);elems(element,2)];       %get nodes of the element
    u_nodes=[u(2*node_element(1,1)-1:2*node_element(1,1));u(2*node_element(2,1)-1:2*node_element(2,1))];    %get displacement at the nodes
    
    u_prim=T(:,:,element)*u_nodes;      %calculate local displacement
    
    epsilon=(u_prim(1,1)-u_prim(3,1))/L(element);   %calculate strain

end

