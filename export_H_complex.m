function export_H_complex(H_complex, filename)
    syms q1 q2 q3 p1 p2 p3
    vars = [q1 q2 q3 p1 p2 p3];
    [coeffs_list, monomials] = coeffs(expand(H_complex), vars);
    fid = fopen(filename, 'w');
    % 每行: e1 e2 e3 e4 e5 e6  re  im
    for i = 1:length(monomials)
        m = monomials(i);
        c = coeffs_list(i);
        exps = zeros(1, 6);
        for k = 1:6
            exps(k) = feval(symengine, 'degree', m, vars(k));
        end
        c_eval = vpa(c, 30);
        re = double(real(c_eval));
        im = double(imag(c_eval));
        fprintf(fid, '%d %d %d %d %d %d %.17g %.17g\n', ...
            exps(1), exps(2), exps(3), exps(4), exps(5), exps(6), re, im);
    end
    fclose(fid);
end