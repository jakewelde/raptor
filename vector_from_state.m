function x = vector_from_state(xs, Rg, th1, th2, xs_d, w, th1d, th2d)
    x = [
        xs;
        reshape(Rg,[9 1]);
        th1;
        th2;
        xs_d;
        w;
        th1d;
        th2d;
    ];
end