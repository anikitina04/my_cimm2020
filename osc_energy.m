function [pot_en, kin_en] = osc_energy(u, v, omega)

    pot_en = 0.5 .* omega^2 .* u.^2;
    kin_en = 0.5 .* v.^2;

end