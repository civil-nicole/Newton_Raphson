function [K_global, K_reduced, T, L] = get_material_stiffness_matrix(nodes, elems, rests)
% This function calculates the material stiffness matrices for the global
% and the reduced system.
%Author: Nicole Widmer
%Last update: 23/10/2021

    %Inputs
    connectivity=elems(:,1:2);                %for every element, indicate the nodes it connects    
    
    A=elems(:,3);                     % vector with sections of elements [mm^2]
    E=elems(:,4);                                % [N/mm^2]
    
    number_elements=size(connectivity,1);       %number of elements
    number_nodes=max(max(connectivity));    %number of nodes
    
    % Unit stiffness matrix of element in local coordinate system
    K_unit=[1 0 -1 0; 
            0 0 0 0; 
           -1 0 1 0; 
            0 0 0 0];
    
    % Geometrical calculations: length and angle of elements
    for element = 1:number_elements                 %loop over every element to calculate the length and the orientation angle of every element
        node1=connectivity(element,1);
        node2=connectivity(element,2);
        coord_node1=nodes(node1,:);     %x and y coordinates of node 1
        coord_node2=nodes(node2,:);     %x and y coordinates of node 2
        
        L(element)= get_length(coord_node1, coord_node2); %[mm] calculation of element length
        theta(element)=get_theta(coord_node1, coord_node2); % calculation of orientation angle of the element
    
        % Element stiffness matrices
        T(:,:,element)=get_transformation_matrix(theta(element));
        
        k=E(element)*A(element)/L(element);           %linear stiffness of the element
        k_element{element}=k*T(:,:,element)'*K_unit*T(:,:,element);      %calculation of the global stiffness matrix of the element        
    end
    
    % Assemble global stiffness matrix
    NDoF=2*number_nodes;
    
    % Assemble all the elements
    [K_global,K_reduced]=assembly_element_matrices(k_element,nodes,elems,rests);

end

