function Ctt = constraint_dtt(revolute, simple, driving, t, q, qt)

r_len = length(revolute);
s_len = length(simple);
d_len = length(driving);

n_constr = 2 * r_len + s_len + d_len;

Ctt = zeros(n_constr, 1);

c_idx = 0;

for r = revolute
    Ctt(c_idx + (1:2), 1) = ...
        revolute_joint_g(r.i, r.j, r.s_i, r.s_j, q, qt);
    c_idx = c_idx + 2;
end