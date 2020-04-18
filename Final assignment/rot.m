function A = rot(phi)
% Rotational transformation matrix for planar case
cp = cos(phi);
sp = sin(phi);
A = [cp, -sp
    sp, cp];