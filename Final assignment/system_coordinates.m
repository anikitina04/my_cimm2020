function q = system_coordinates(bodies)

b_len = length(bodies);

q = zeros(3 * b_len, 1);

for ii = 1 : b_len
    q(body_idx(ii)) = bodies(ii).q;
end