function F = evolut(ti, x, r, tau)
syms t real;
sigmx = [ 0 1 ; 
          1 0];
sigmz = [ 1 0 ; 
         0 -1];
sigmy = [0 -1i; 
         1i 0];
q = [ 1 0 ; 
      0 1];
cnot = [1 0 0 0; 
        0 1 0 0; 
        0 0 0 1; 
        0 0 1 0];
H_s = [1i*r 1; 
       1 -1i*r] ;
% read supplementary of the article https://science.sciencemag.org/content/364/6443/878  
M = findM(t, r, tau);

%determine \eta as eq.15 in supplementary https://science.sciencemag.org/content/364/6443/878
N_pre = sqrtm(M - eye(2));

N = [piecewise(t == 0,sqrt(subs(M(1),t,0)- 1), N_pre(1)) piecewise(t == 0,0, N_pre(3)); 
    piecewise(t == 0,0, N_pre(2)) piecewise(t == 0,sqrt(subs(M(4),t,0)- 1), N_pre(4))];
% d(\eta)/dt
dif_N0 = [limit((N(1) - subs(N(1),t,0))/t,t,0,'right') limit((N(3) - subs(N(3),t,0))/t,t,0,'right'); 
 limit((N(2) - subs(N(2),t,0))/t,t,0,'right') limit((N(4) - subs(N(4),t,0))/t,t,0,'right')];

dif_N_pre = diff(N);

diff_N = [piecewise(t == 0, dif_N0(1), dif_N_pre(1)) , piecewise(t == 0, dif_N0(3), dif_N_pre(3)) ;
 piecewise(t == 0, dif_N0(2), dif_N_pre(2)) piecewise(t == 0, dif_N0(4), dif_N_pre(4))];

%M^(-1)
inver_M = eye(2)/M;


% eq. 21
J = (H_s + (1i*diff_N + N*H_s)*N)*inver_M;
W = 1i*(H_s*N-N*H_s-1i*diff_N)*inver_M;

%eq. 20
H_sa = kron(J,eye(2))+kron(W,sigmz);
H_sai = subs(H_sa,t,ti);
F = double(-1i*H_sai*x);

end

