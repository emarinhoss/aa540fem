function [kxx,kxy,kyx,kyy,f]=conductivity_and_forcing(x,y)

%%%
% conductivity tensor definition and
% forcing vector analitical equations
% k(x,y)
% f(x,y)
%
% NOTE: use lower case for 'x' and 'y'
%%%

kxx = 1;
kxy = 0;
kyx = 0;
kyy = 1;

f=0;
%f = (sin(pi*x)*sin(pi*y))^(40);