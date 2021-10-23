function [T]=get_transformation_matrix(theta)
%This function calculates the transformation matrix for an element based on
%the orientation (i.e. angle theta)
%Author: Nicole Widmer
%Last update: 23/10/2021

    c=cosd(theta);
    s=sind(theta);
    
    T = [c s 0 0; 
        -s c 0 0; 
        0 0 c s; 
        0 0 -s c];

end