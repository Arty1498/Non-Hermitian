function F = mineigval1(t, r)
H_s = [1i*r 1; 
       1 -1i*r];
% R is a matrix from eq. (12) in the supplementary https://science.sciencemag.org/content/364/6443/878
R = 2*expm(-1i*t*H_s')*expm(1i*t*H_s);
v = eig(R);
%return the first eigenvalue
F = v(1);
end
