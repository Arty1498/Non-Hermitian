function F = findM(t, r, tau)
%find the minimum eigenvalue on the time interval [0, tau]
[x1,fval1] = fminbnd(@(x) mineigval1(x, r), 0, tau);
[x2,fval2] = fminbnd(@(x) mineigval2(x, r), 0, tau);
myu = min(fval1, fval2);
myu = 0.99*myu;


H_s = [1i*r 1; 
      1 -1i*r];
%determine the matrix M
F = 2*expm(-1i*t*H_s')*expm(1i*t*H_s)/myu;
end

