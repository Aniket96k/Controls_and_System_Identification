function [K,Gzw] = Hinf(A,B1,B2,C1,D11,D12,gamma)
%HINF Summary of this function goes here
%   Detailed explanation goes here

%Sizing
n = size(A,1) ;
m1 = size(B1,2) ;
m2 = size(B2,2) ;
p1 = size(C1,1) ;
Im = eye(m1) ;
Ip = eye(p1) ;
In = eye(n) ;

%Set LMI
setlmis([]) 

%Matrix Variables
% gamma = lmivar(1,[1 1]) ;
X = lmivar(1,[n 1]) ;
W = lmivar(2, [m2 n]) ;

%LMI Terms
lmiterm([1 1 1 X],A,1,'s');    %X'A'   (1,1)
lmiterm([1 1 1 W],B2,1,'s');   %W'B2'  (1,1)
lmiterm([1 1 2 0],B1);       %B1     (1,2)
lmiterm([1 1 3 X],1,C1');   %X'C1'  (1,3)
lmiterm([1 1 3 -W],1,D12');  %W'D12' (1,3)
lmiterm([-1 2 2 0],gamma); %-yI    (2,2)
lmiterm([1 2 3 0],D11');     %D11'   (2,3)
lmiterm([-1 3 3 0],gamma); %-yI    (3,3)
lmiterm([-2 1 1 X],1,1);     %X>0

lmi_sys = getlmis ;

%Feasibility Problem
[tmin,xfeas] = feasp(lmi_sys) ;
X = dec2mat(lmi_sys,xfeas,X) ;
W = dec2mat(lmi_sys,xfeas,W) ;

K = W*inv(X) ;
sys = ss(A+B2*K,B1,C1+D12*K,D11) ;
Gzw = hinfnorm(sys) ;

end

