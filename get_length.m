function [L] = get_length(coordinates_node1, coordinates_node2)


L=sqrt((coordinates_node1(1,1)-coordinates_node2(1,1))^2+(coordinates_node1(1,2)-coordinates_node2(1,2))^2);

end