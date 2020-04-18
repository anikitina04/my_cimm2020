function C = constraint(revolute, simple, driving, t, q)

r_len = length(revolute);
s_len = length(simple);
d_len = length(driving);

n_constr = 2 * r_len + s_len + d_len;

C = zeros(n_constr, 1);

c_idx = 0;
for r = revolute
    C(c_idx + (1:2)) = revolute_joint(r.i, r.j, r.s_i, r.s_j, q);
    c_idx = c_idx + 2;
end

for s = simple
    c_idx = c_idx + 1;
    C(c_idx) = simple_joint(s.i, s.k, s.c_k, q);
end

for d = driving
    c_idx = c_idx + 1;
    C(c_idx) = driving_joint(d.i, d.k, d.d_k, t, q);
end