function x = vector_from_state(xs, Rq, Rg, xs_d, Om, w)
    x = [
        xs;
        reshape(Rq,[9 1]);
        reshape(Rg,[9 1]);
        xs_d;
        Om;
        w;
    ];
end