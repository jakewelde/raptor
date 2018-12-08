function draw_robot(state)
    global mq_ mg_ Lg_ Le_ e1
    global compute_Rq_state
    
    [xs, Rg, th1, th2, xs_d, w, th1d, th2d] = state_from_vector(state);
    
    Rq = compute_Rq_state(Rg,th1,th2);
    xg = xs+Lg_*mq_/(mg_+mq_)*Rg*e1;
    xq = xg - Lg_*Rg*e1;
        
    
    
    coords = .5*[
        1 -1  0  0
        0  0  1 -1
        0  0  0  0
      ];
    hold on
    for row = coords
        R = .3 ;    % Radius of circle 
        teta=0:0.01:2*pi ;
        xyz = Rq*[
            row(1) + R*cos(teta);
            row(2) + R*sin(teta) ;
            row(3) + zeros(size(teta));
        ];
        patch(xq(1) + xyz(1,:),xq(2) + xyz(2,:),xq(3) + xyz(3,:),[.7 .8 .7])
    end
    alpha(.4)

    bar1 = xq*[1 1]+Rq*coords(:,1:2); % quad frame x
    bar2 = xq*[1 1]+Rq*coords(:,3:4); % quad frame y
    
    xe = xq+(xg-xq)*Le_/Lg_;
    
    bar3 = [xq xe];                   % arm
    bar4 = xe*[1 1]+.1*Rg*[           % cross piece of gripper
        0  0;
        1 -1;
        0  0;
    ];
    bar5 = xe*[1 1]+.1*Rg*[           % left finger
        0  2;
        -1 -1;
        0  0;
    ];              
    bar6 = xe*[1 1]+.1*Rg*[           % right finger
        0  2;
        1 1;
        0  0;
    ];

    plot3(bar1(1,:),bar1(2,:),bar1(3,:),'color',[.2 .4 .2],'LineWidth',2)
    plot3(bar2(1,:),bar2(2,:),bar2(3,:),'color',[.2 .4 .2],'LineWidth',2)
    plot3(bar3(1,:),bar3(2,:),bar3(3,:),'color',[.2 .2 .4],'LineWidth',2)
    plot3(bar4(1,:),bar4(2,:),bar4(3,:),'color',[.2 .2 .4],'LineWidth',2)
    plot3(bar5(1,:),bar5(2,:),bar5(3,:),'color',[.2 .2 .4],'LineWidth',2)
    plot3(bar6(1,:),bar6(2,:),bar6(3,:),'color',[.2 .2 .4],'LineWidth',2)

    
    
    quiver3(xq(1),xq(2),xq(3),Rq(1,1),Rq(2,1),Rq(3,1),'AutoScale','off','color',[1 0 0],'LineWidth',1.5,'MaxHeadSize',.25);
    quiver3(xq(1),xq(2),xq(3),Rq(1,2),Rq(2,2),Rq(3,2),'AutoScale','off','color',[0 1 0],'LineWidth',1.5,'MaxHeadSize',.25);
    quiver3(xq(1),xq(2),xq(3),Rq(1,3),Rq(2,3),Rq(3,3),'AutoScale','off','color',[0 0 1],'LineWidth',1.5,'MaxHeadSize',.25);

    Rg = 1/2*Rg; % visualization scaling
    quiver3(xg(1),xg(2),xg(3),Rg(1,1),Rg(2,1),Rg(3,1),'AutoScale','off','color',[1 0 0],'LineWidth',1.5,'MaxHeadSize',.5);
    quiver3(xg(1),xg(2),xg(3),Rg(1,2),Rg(2,2),Rg(3,2),'AutoScale','off','color',[0 1 0],'LineWidth',1.5,'MaxHeadSize',.5);
    quiver3(xg(1),xg(2),xg(3),Rg(1,3),Rg(2,3),Rg(3,3),'AutoScale','off','color',[0 0 1],'LineWidth',1.5,'MaxHeadSize',.5);

    plot3(xs(1),xs(2),xs(3),'o',...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[.8 .8 1],...
        'MarkerSize',7);
end