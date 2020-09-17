function F = mineigval2(t, r)
H_s = [1i*r 1; 1 -1i*r];
R = 2*expm(-1i*t*H_s')*expm(1i*t*H_s);
v = eig(R);
F = v(2);
end