function [xs, Rg, th1, th2, xs_d, w, th1d, th2d] = state_from_vector(x)
    xs = x(1:3);
    Rg = reshape(x(4:4+8),[3,3]);
    th1 = x(13);
    th2 = x(14);
    xs_d = x(15:17);
    w = x(18:20);
    th1d = x(21);
    th2d = x(22);
end