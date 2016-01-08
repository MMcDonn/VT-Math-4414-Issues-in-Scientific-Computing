function int_val = quadrature_trap_comp_f_vec_nonunif(f_vec, t)
    int_val = 0;
for i = 1:length(t)-1
    int_val = int_val + ((t(i+1) - t(i))/2)*(f_vec(i) + f_vec(i+1));
end