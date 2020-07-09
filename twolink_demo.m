% manipulability and force ellipsoid of two links
%
% Author: Yichi Iino
%
% Refrerence:
%
%

clear;

L1 = 0.4; L2 = 0.4;

for iPos = 1:9
    
    % set link relative angle 
    theta1 = pi/2*iPos/10;
    theta2 = 2*(pi/2 - theta1);
    
    % set link position
    x1 = L1*sin(theta1); y1 = L1*cos(theta1);
    x2 = L1*sin(theta1) + L2*sin(theta1 + theta2);
    y2 = L1*cos(theta1) + L2*cos(theta1 + theta2);
    
    % plot link
    figure(1)
    subplot(2,1,1)
    plot([0, x1], [0, y1], 'k'); hold on;
    plot([x1, x2], [y1, y2], 'k');
    
    % set the endpoint coordinate
    x = L1*sin(theta1) + L2*sin(theta1 + theta2);
    y = L1*cos(theta1) + L2*cos(theta1 + theta2);
    
    % set Jacobian matrix
    J = [L1*cos(theta1) + L2*cos(theta1 + theta2) L2*cos(theta1 + theta2); ...
        -L1*sin(theta1) - L2*sin(theta1 + theta2) -L2*sin(theta1 + theta2)];
    
    J1 = inv(J)'*inv(J);
    J2 = J*J.';
    
    % determine eigen values and vectors
    [V1, D1] = eig(J1);
    [V2, D2] = eig(J2);
    
    % plot manipulability ellipsoid
    t = 0:0.05*pi:2*pi;
    a = 1/sqrt(D1(1,1)); b = 1/sqrt(D1(2,2));
    
    xt = a*cos(t).*V1(:,1) + b*sin(t).*V1(:,2) + [x; y];
    plot(xt(1, :), xt(2, :))
    
    % plot links
    subplot(2,1,2)
    plot([0, x1], [0, y1], 'k'); hold on;
    plot([x1, x2], [y1, y2], 'k');
    
    % plot force ellipse
    t = 0:0.05*pi:2*pi;
    a2 = 1/sqrt(D2(1,1)); b2 = 1/sqrt(D2(2,2));
    a2 = a2*0.05; b2 = b2*0.05;
    
    xt2 = a2*cos(t).*V2(:,1) + b2*sin(t).*V2(:,2) + [x; y];
    plot(xt2(1,:), xt2(2, :))
    
end
    
    
    