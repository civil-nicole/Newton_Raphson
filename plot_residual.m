% Plot residual force
% This script plots the evolution of the residual force for all the
% iterations
%Author: Nicole Widmer
%Last update: 23/10/2021

%read input file
R=readmatrix("R.txt");

%plotting
figure()
hold on
plot(R(:,1)'/1000)
plot(R(:,2)'/1000)
xlabel('Iteration')
ylabel('Residual force R [kN]')
legend(['x-direction';'y-direction'],'Location', 'Northeast')
title('Residual force')
hold off