function [K_global,K_red] = assembly_element_matrices(K_bars,nodes,elems,rests)
% This function assembles the stiffness matrices of the different elements
% into the global stiffness matrix. It returns the global stiffness matrix
% and the stiffness matrix of the reduced system.
%Author: Nicole Widmer
%Last update: 23/10/2021

    number_bars=size(elems,1);
    number_nodes=size(nodes,1);

    connectivity=elems(:,1:2);                %for every bar, indicate the nodes it connects
    i=1;
    for node=1:size(nodes,1)
       for j=1:2
           degrees_of_freedom(node,j)=rests(i);
           i=i+1;
       end
    end

    
    % assemble global stiffness matrix
    NDoF = 2*size(nodes, 1);
    K_global=zeros(NDoF);

    for m=1:number_bars 
        % Node 1 of element
        N1=elems(m,1);
        % Node 2 of element
        N2=elems(m,2);
        % find indices of DOFs
        i1=(N1*2-1); i2=N1*2; % indeces of DOFs of Node 1
        i3=(N2*2-1); i4=N2*2; % indeces of DOFs of Node 2
        ii=[i1 i2 i3 i4];

        % Expand element stiffness matrix in global CS by the missing DOFs
        Ki = zeros(size(K_global));
        for k=1:4
            for l=1:4
                Ki(ii(k),ii(l))=K_bars{m}(k,l);
            end
        end
        K_global=K_global+Ki;
        
    end



    % reduced system
    counter=1;

    for node=1:number_nodes
       %x-direction
       if degrees_of_freedom(node,1)==1     %if DOF is blocked, set displacement to 0
            u(2*node-1,1)=0;
       else
           u(2*node-1,1)=1;                 %if DOF is free, set displacement to 1
           DOF_line(counter,1)=2*node-1;    %store which line in global stiffness matrix corresponds to nodes and direction of reduced system
           counter=counter+1;
       end

       %y-direction
       if degrees_of_freedom(node,2)==1
            u(2*node,1)=0;                  %if DOF is blocked, set displacement to 0 
       else
           u(2*node,1)=1;                   %if DOF is free, set displacement to 1
           DOF_line(counter,1)=2*node;      %store which line in global stiffness matrix corresponds to nodes and direction of reduced system
           counter=counter+1;
       end

    end
    
    counter_R=1;
    for i=1:size(u,1)

        if u(i,1)==1
            nodes_red(counter_R,1)=(i+1)/2;         %array with all the nodes in reduced system (y-direction for node 1: node = 1.5)         %vector with reactions for reduced system
            counter_R=counter_R+1;
        end

    end
    
    
    %stiffness matrix for reduced system
    counter_K=1;
    for i=1:size(nodes_red,1)
        node=nodes_red(i,1);

        %store the columns for reduced system
        if round(node)==1           %integer node = x-direction
            K_red_interm(:,counter_K)=K_global(:,round(2*node));        
        else
            K_red_interm(:,counter_K)=K_global(:,round(2*node-1));
        end
        counter_K=counter_K+1;
    end
    
    % consider only the rows corresponding to free DOFs
    for i=1:size(DOF_line,1)
            K_red(i,:)=K_red_interm(DOF_line(i),:);
    end

end