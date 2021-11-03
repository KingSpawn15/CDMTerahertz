% Based on the paper: 'Electro-optic transceivers for terahertz-wave
% applications' Q. Chen, M. Tani, Zhiping Jiang, and X.-C. Zhang JOSA B 2001
% pages 826-827

clc; close all; clear all;

%% Lattice rotation matrix
% [x',y',z'] - Lab ref. frame. x' is the <111> direction (normal to surface)
% [x,y,z] - Cryatal ref. frame. x is the <100> direction etc.

xprime = [1;1;1];
xprime = xprime./sqrt(sum(abs(xprime).^2));

zprime = [-1;-1;2];
zprime = zprime./sqrt(sum(abs(zprime).^2));

yprime = cross(zprime,xprime);
yprime = yprime./sqrt(sum(abs(yprime).^2));

%From lab ref. frame (x'y'z') to crystal ref. frame (xyz)
RLab2Crystal = [xprime yprime zprime]; 

%%
syms theta

%In lab ref. frame (x'y'z'). Assuming AOI=0 and linear polarization in
%y'-z' plane
% theta is the angle from z' axis. CCW is poisitive
Eprime = [0;-sin(theta);cos(theta)];

%In crystal ref. frame (xyz)
E = simplify(RLab2Crystal*Eprime);
P = simplify(2.*[E(2)*E(3);E(3)*E(1);E(1)*E(2)]);

Pprime = simplify(inv(RLab2Crystal)*P);