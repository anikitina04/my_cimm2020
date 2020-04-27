%% Second derivative of revolute joint

function g = revolute_joint_g(i, j, s_i, s_j, q, qd)

idx_i = body_idx(i);
phi_i = q(idx_i(3));
phi_i_d = qd(idx_i(3));

idx_j = body_idx(j);
phi_j = q(idx_j(3));
phi_j_d = qd(idx_j(3));

g = (rot(phi_i) * s_i * phi_i_d^2) - (rot(phi_j) * s_j * phi_j_d^2);
end