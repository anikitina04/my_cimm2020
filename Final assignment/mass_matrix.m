function M = mass_matrix(bodies)

b_len = length(bodies);

M = zeros(b_len * 3);

for ii = 1:b_len
    b = bodies(ii);
    idx = body_idx(ii);
    M(idx, idx) = diag([b.m, b.m, b.Ic]);
end