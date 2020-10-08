function F = mineigval1(t, r)
H_s = [1i*r 1; 
       1 -1i*r];
% M(t) eq. (12)
M = 2*expm(-1i*t*H_s')*expm(1i*t*H_s);
v = eig(M);
%return the first eigenvalue
F = v(1);
end
