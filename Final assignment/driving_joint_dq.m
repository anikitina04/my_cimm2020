function Cq = driving_joint_dq(k)
% i - body id
% k = 1, 2, 3 for x, y and phi
% d_k - function of time that define coordinate value
% q - coordinate vector
Cq = zeros(1, 3);
Cq(k) = 1;