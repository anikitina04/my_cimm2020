% Script example 1
r_3 = [1.2; 2.5];
phi_3 = 5.67;
sB_local = [-1.8; 1.3];

% global components of vector sB3
sB_3 = rot(phi_3) * sB_local;
disp(sB_3);

rB = r_3 + sB_3;
disp(rB);