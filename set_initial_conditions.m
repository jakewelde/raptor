% Hover
Rq0 = eye(3);

% Mix
% Rg0 = axisangle(e1,pi/8)*axisangle(e2,.45*pi);
% T2 direction
Rg0 = axisangle(e2,pi/8);

x0 = vector_from_state(...
    [0;0;0],Rq0,Rg0,...
    [3;0;0],[0;0;0],Rg0.'*[0;0;0]...
);


% T1 direction
% Rg0 = axisangle(e1,pi/6)*axisangle(e2,pi/2);

% Down
% Rg0 = axisangle(e2,pi/2);




% twist
% x0 = vector_from_state(...
%     [0;0;0],Rq0,Rg0,...
%     zeros(3,1),[0;0;.1],Rg0.'*[0;0;.1]...
% );
