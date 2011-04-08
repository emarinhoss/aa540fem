function [xi,wi] = gauss_legendre_quad(o)

switch o
    case 2
        xi = [-(sqrt(3))^-1 (sqrt(3))^-1];
        wi = [1 1];
    case 3
        xi = [-sqrt(15)/5 0 sqrt(15)/5];
        wi = 1/9*[5 8 5];
    case 4
        xi = [-sqrt((3+2*sqrt(6/5))/7) -sqrt((3-2*sqrt(6/5))/7) ...
               sqrt((3-2*sqrt(6/5))/7) sqrt((3+2*sqrt(6/5))/7)];
        wi = [(18-sqrt(30))/36 (18+sqrt(30))/36 ...
              (18+sqrt(30))/36 (18-sqrt(30))/36];
    otherwise
        xi = 0;
        wi = 2;
end