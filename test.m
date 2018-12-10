        clc

        velocity = @(C) C(17:24)*(0:7).' + e3.'*Ls_ * compute_Rg_angles(sum(C(25:32)),sum(C(33:40)),sum(C(41:48)))*...
        hat(compute_w_angles(sum(C(1:8)),C(25:32)*(0:7).',sum(C(9:16)),C(33:40)*(0:7).',C(41:48)*(0:7).'))...
        *e1 - zdf(3),
    
    zdf = 0*ones(3,1);
    C = [
        zeros(8*3,1);
        zeros(8,1);
        [1 1 0 0 0 0 0 0].';
        zeros(8,1);
    ].'

velocity(C)