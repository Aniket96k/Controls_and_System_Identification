% clc, clear all, close all, format compact

%% a
A = A_cap ;
B2 = B_cap ;
B1 = [0 0 1 1]' ;
C1 = eye(4) ;
A_rows = size(A,1) ;
B1_cols = size(B1,2) ;
B2_cols = size(B2,2) ;
C_rows = size(C1,1) ;
D11 = zeros(C_rows,B1_cols) ;
D12 = zeros(C_rows,B2_cols) ;
gamma = .5 ;

[K,Gzw] = Hinf(A,B1,B2,C1,D11,D12,gamma)

%% b
sys = ss(A+B2*K,B1,C1+D12*K,D11) ;

[y,t,x] = impulse(sys,5) ;
figure(1)
impulse(sys)

%% c
figure(2)
sigma(sys)
fprintf('The system attenuates the gain as seen by the negative DB (<1 abs) gain values.\n')

%% d
u = K*x' ;
L2_u1 = norm(u(1,:)) 
L2_u2 = norm(u(2,:))
fprintf('The first control signal requires more energy.\n')

figure(3)
plot(t,(K*x')')
title('Control Signals'); xlabel('time (s)'); ylabel('u(t)')
legend('u1','u2')