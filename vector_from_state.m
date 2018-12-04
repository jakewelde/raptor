function x = vector_from_state(xs, Rq, Rg, xs_d, Om, w, xb, xb_d)
    x = [
        xs;
        reshape(Rq,[9 1]);
        reshape(Rg,[9 1]);
        xs_d;
        Om;
        w;
        xb;
        xb_d;
    ];
end