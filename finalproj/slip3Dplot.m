%% SLIP3D Visualization Script
clear; clc; close all;

% Load saved data
load('SLIP3D_data.mat');

% traj: 3×N (rows = [h0; vx; vy])
% K: gains
% X0_star, u0_star: optimal control pair for desired gait timing

%% Plot apex height and velocities
figure;
subplot(3,1,1);
plot(traj(1,:), 'o-', 'LineWidth', 1.5);
xlabel('Step #'); ylabel('Apex Height (m)');
title('Apex Height per Step'); grid on;

subplot(3,1,2);
plot(traj(2,:), 'o-', 'LineWidth', 1.5);
xlabel('Step #'); ylabel('Forward Velocity (m/s)');
title('Forward Velocity per Step'); grid on;

subplot(3,1,3);
plot(traj(3,:), 'o-', 'LineWidth', 1.5);
xlabel('Step #'); ylabel('Lateral Velocity (m/s)');
title('Lateral Velocity per Step'); grid on;

sgtitle('3D SLIP Apex States Across Steps');


%% Print summary
if exist('u0_star','var') && exist('X0_star','var')
    fprintf('Optimized Periodic Gait:\n');
    fprintf('  Apex height h0 = %.3f m\n', X0_star(1));
    fprintf('  vx = %.3f m/s, vy = %.3f m/s\n', X0_star(2), X0_star(3));
    fprintf('  θ = %.2f°, φ = %.2f°, ks1 = %.1f, ks2 = %.1f\n', ...
        rad2deg(u0_star(1)), rad2deg(u0_star(2)), u0_star(3), u0_star(4));
end