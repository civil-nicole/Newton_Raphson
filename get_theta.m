function [theta] = get_theta(coordinates_node1, coordinates_node2)

theta=atand((coordinates_node2(1,2)-coordinates_node1(1,2))/(coordinates_node2(1,1)-coordinates_node1(1,1)));

end