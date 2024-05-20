k = 1 ;
Ts = .002 ;
T = 15;
Tt = 5;
wr = 2*pi/T ;
N = floor(100/wr) ;
w0 = 7/8*2*pi ;
w1 = 10*w0;
tf = 's';
% filter = s/(s/w1+1)^2;
[A1,B1,C1,D1] = tf2ss([1 0], [1/w1^2 sqrt(2)/w1 1]);
wN = pi/Ts ;
t = 0:Ts:T-Ts ;
w = (0:wr:2*wN-wr)' ;
Beq = .004;
g = 9.81;
Jeq = .0035842;
Jm = 3.87E-7;
Kg = 70;
Km = .00767;
Kt = Km;
L = .1675;
m = .125;
r = .215;
Rm = 2.6;
zeta_g = .9;
zeta_m = .69;
a = Jeq+m*r^2+zeta_g*Kg^2*Jm;
b = m*L*r;
c = 4/3*m*L^2;
d = m*g*L;
e = Beq+(zeta_m*zeta_g*Kt*Kg^2*Km)/Rm;
f = (zeta_m*zeta_g*Kt*Kg)/Rm;
M = [a -b; -b c];
G = [0 0; 0 d];
Di = [e 0; 0 0];
Bu = [f; 0];

MiG = M\G;
MiD = M\Di;
MiB = M\Bu;
A = [0 0 1 0; 0 0 0 1; MiG MiD];
B = [0; 0; MiB];
C = eye(4);
D = [0 0 0 0]';
R = 1;
x1_nom = 1;
x2_nom = 100;
x3_nom = 100;
x4_nom = 100;
Q11 = 1/(x1_nom)^2;
Q22 = 1/(x2_nom)^2;
Q33 = 1/(x3_nom)^2;
Q44 = 1/(x4_nom)^2;
Q = diag([Q11; Q22; Q33; Q44]);