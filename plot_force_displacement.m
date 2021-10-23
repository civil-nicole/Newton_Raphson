% Plot force-displacement response
%This script plots the force-displacement response (from U.txt and F.txt) for all degrees of freedom
%for a specified node
%Author: Nicole Widmer
%Last update: 23/10/2021

clc
clear all

%define note to be plotted
node_to_plot=2;

% read input files
U=readmatrix("U.txt");
F=readmatrix("F.txt");


% degrees of freedom corresponding to node to be plotted
DOF_to_plot=[2*node_to_plot-1,2*node_to_plot];

%plot force displacement curves for all degrees of freedom
for i=1:size(DOF_to_plot,2)
    figure(i)
    hold on
    plot([0 U(DOF_to_plot(i),:)],[0 F(DOF_to_plot(i),:)/1000])
    set(gca,'Ydir', 'reverse')
    set(gca,'Xdir', 'reverse')
    
    if mod(DOF_to_plot(i),2)~=0
    xlabel('Nodal displacement U1 [m]')
    ylabel('Nodal force P1 [kN]')
    else
    xlabel('Nodal displacement U2 [m]')
    ylabel('Nodal force P2 [kN]')
    end

    title(['Newton-Raphson: Force-displacement at Node ' num2str(node_to_plot)])
    hold off
end

