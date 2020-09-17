function F = root(x, U)
cnot = [1 0 0 0; 
       0 1 0 0; 
       0 0 0 1; 
       0 0 1 0];
%General form of the U3 unitary matrix
c1 = [cos(x(1)/2) -exp(x(9)*1i)*sin(x(1)/2)  ; 
    exp(x(17)*1i)*sin(x(1)/2) exp((x(9)+x(17))*1i)*cos(x(1)/2)];
c2 = [cos(x(2)/2) -exp(x(10)*1i)*sin(x(2)/2) ; exp(x(18)*1i)*sin(x(2)/2) exp((x(10)+x(18))*1i)*cos(x(2)/2)];
c3 = [cos(x(3)/2) -exp(x(11)*1i)*sin(x(3)/2) ; exp(x(19)*1i)*sin(x(3)/2) exp((x(11)+x(19))*1i)*cos(x(3)/2)];
c4 = [cos(x(4)/2) -exp(x(12)*1i)*sin(x(4)/2) ; exp(x(20)*1i)*sin(x(4)/2) exp((x(12)+x(20))*1i)*cos(x(4)/2)];
c5 = [cos(x(5)/2) -exp(x(13)*1i)*sin(x(5)/2) ; exp(x(21)*1i)*sin(x(5)/2) exp((x(13)+x(21))*1i)*cos(x(5)/2)];
c6 = [cos(x(6)/2) -exp(x(14)*1i)*sin(x(6)/2) ; exp(x(22)*1i)*sin(x(6)/2) exp((x(14)+x(22))*1i)*cos(x(6)/2)];
c7 = [cos(x(7)/2) -exp(x(15)*1i)*sin(x(7)/2) ; exp(x(23)*1i)*sin(x(7)/2) exp((x(15)+x(23))*1i)*cos(x(7)/2)];
c8 = [cos(x(8)/2) -exp(x(16)*1i)*sin(x(8)/2) ; exp(x(24)*1i)*sin(x(8)/2) exp((x(16)+x(24))*1i)*cos(x(8)/2)];

%our cicuit 
C = kron(c7,c8)*cnot*kron(c5,c6)*cnot*kron(c3,c4)*cnot*kron(c1,c2);

%a function of the accuracy of the approximation
T = C*U' - eye(4);

F = [real(T(1));
    imag(T(1));
    real(T(2));
    imag(T(2));
    real(T(3));
    imag(T(3));
    real(T(4));
    imag(T(4));
    real(T(5));
    real(T(6));
    imag(T(6));
    real(T(7));
    imag(T(7));
    real(T(8));
    imag(T(8));
    imag(T(9));
    real(T(10));
    real(T(11));
    imag(T(11));
    real(T(12));
    imag(T(12));
    imag(T(14));
    real(T(16));
    imag(T(16))];
end





