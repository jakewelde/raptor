function [xs, Rq, Rg, xs_d, Om, w] = state_from_vector(x)
    xs = x(1:3);
    Rq = reshape(x(4:4+8),[3,3]);
    Rg = reshape(x(13:13+8),[3,3]);
    xs_d = x(22:24);
    Om = x(25:27);
    w = x(28:30);
end