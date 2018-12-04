while true
    cx = rand([8,1])-.5; cy = rand([8,1])-.5; cz = rand([8,1])-.5;

    time_range = 0:.01:1;
    positions = zeros(3,size(time_range,2));
    for i=1:size(time_range,2)
        positions(:,i) = [
            [1 0 0 0 0]*compute_derivatives(cx,time_range(i),1);
            [1 0 0 0 0]*compute_derivatives(cy,time_range(i),1);
            [1 0 0 0 0]*compute_derivatives(cz,time_range(i),1)
        ];
    end

    plot3(positions(1,:),positions(2,:),positions(3,:))
    axis equal;
    pause;
    
end