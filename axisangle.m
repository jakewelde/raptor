function R = axisangle(u,a)

    u = u / norm(u);
    ux = u(1);
    uy = u(2);
    uz = u(3);
    
    ca = cos(a);
    sa = sin(a);
    
    R = [
        ca + ux^2*(1-ca)     ux*uy*(1-ca)-uz*sa     ux*uz*(1-ca)+uy*sa;
        uy*ux*(1-ca)+uz*sa   ca+uy^2*(1-ca)         uy*uz*(1-ca)-ux*sa;
        uy*uz*(1-ca)-uy*sa   uz*uy*(1-ca)+ux*sa     ca+uz^2*(1-ca)
    ];

end