function Cq = simple_joint_dq(k)
% i - body id
% k = 1, 2, 3 for x, y and phi
% c_k - value of coordinate to keep all the time
% q - coordinate vector

Cq = zeros(1, 3);
Cq(k) = 1;