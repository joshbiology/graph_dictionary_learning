function OMP_wrapper(dict_path, signal_path, out_name, t_sparsity)
%=============================================
% Wrapper for OMP
% Convenience function for Pan et al., 2021


rng('default')



dict = readmatrix(dict_path);

dict = normcols(dict);
signal = readmatrix(signal_path);


g_mat = dict'*dict;


out = omp(dict, signal, g_mat, t_sparsity);

writematrix(out, out_name);
