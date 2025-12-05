%% 11/27 gait library for updated bdx droid values no LCM
% copy of gaitlib-11-26.m
% has new slip_return_map version w/ 4 event detections
% WIP: updating slip_return_map for position/time continuity between steps? 

clear; clc;
run ../setup.m 

%%
control_freq = 500; % control frequency in Hz
rate_ctrl = rateControl(control_freq);
dt = 1 / control_freq;
steps = 1000; % Tot steps for the simulation

% %% SLIP Parameters
params.M = 2.107141;               % effective point mass
params.g = [0; 0; -9.81];
params.l0 = 0.2019;           % rest spring leg length, at TD l0 = lh
params.lh = 0.2019;            % humanoid virtual leg length used to map to SLIP leg
params.yhip = 0.035;        % zero for testing. lateral hip offset (left hip at y=0.035, right hip at y=-0.035)
params.th0 = deg2rad(25);   % init TD angle guess
params.ks0 = 1.5;             % init stiffness guess (kN/m)
params.tf = 2.0;            % single step time interval
% sim parameters
N = 5; % number of steps
params.X0 = [0.21; 0.6; 0]; % starting X for sim loop

% % robot physical parameters (for dynamics computation)
% params.l = [0.2; 0.0; 0.4; 0.4; 0.03];  % [l_hip_roll; l_hip_pitch; l_thigh; l_shin; l_foot]
% params.M_masses = [1.685713; 0.158036; 0.025867; 0.026811];  % [M_trunk; M_thigh; M_shin; M_foot]
% params.I_inertias = [0.01; 0.005; 0.005; 7.0e-05];  % [I_trunk; I_thigh; I_shin; I_foot]
% params.params = [params.g(3); params.l; params.M_masses; params.I_inertias];

%% build gait library 
vx_range = linspace(0.5, 2.5, 31); % range of desired forward velocities
X0_stars = zeros(3, length(vx_range));
u0_stars = zeros(4, length(vx_range));
K_all = cell(1,length(vx_range));
fprintf('Generating gait library...\n');
for i = 1:length(vx_range)
    vx_des = vx_range(i);
    % speed-dependent stiffness guess: ks0 = a*vx_des + b
    % for vx = 0.5 m/s: ks0 = 1.5 kN/m
    % for vx = 2.5 m/s: ks0 = 5 kN/m
    params.ks0 = 5.0 + (5.0 - 1.5) * (vx_des - 0.5) / (2.5 - 0.5);
    % speed-dependent TD angle guess: th0 should increase with speed
    % for range 0.5-2.5 m/s, use linear interpolation: 22° to 24°
    params.th0 = deg2rad(22.0 + (24.0 - 22.0) * (vx_des - 0.5) / (2.5 - 0.5));
    X0 = [0.25; vx_des; 0]; % initial guess apex state [h, vx, vy], h and vy are guesses to be adjusted, vx is desired vel for entire traj
    
    fprintf('\n--- Speed = %.2f m/s ---\n', vx_des);
    [X0_star, u0_star] = find_periodic_gait(X0, params);
    fprintf('Periodic gait found for %.2f m/s:\n', vx_des);
    fprintf('  Apex height h0     = %.4f m\n', X0_star(1));
    fprintf('  θ = %.2f°, ks = %.3f kN/m\n', rad2deg(u0_star(1)), u0_star(3));
    fprintf('  Lateral velocity vy = %.4f m/s\n', X0_star(3));
    fprintf('  φ     = %.3f deg\n', rad2deg(u0_star(2)));
    fprintf('--------------------------------------------------\n');

    K = compute_deadbeat(X0_star, u0_star, params)
    
    X0_stars(:,i) = X0_star;
    u0_stars(:,i) = u0_star;
    K_all{i} = K;
end
save('SLIP3D_BDX_gait_library.mat','vx_range','X0_stars','u0_stars','K_all','params');
fprintf('\nGait library generation complete.\nSaved to SLIP3D_BDX_gait_library.mat\n');

%% Simulation test
slip_states = zeros(3, N); % to see slip states at each step.
X0 = params.X0; 
pos_log = [0; 0; X0(1)];   % (3x?) for position plot. start at [0, 0, h0]
t_log = [0];               % (1x?) time for position plot.
fprintf('Starting simulation...\n');
for n = 1:N
    % retrieve closest K corresponding to vx from "library":
    fprintf('Retrieving from library...\n');
    fprintf('Step %d: X0 = [%.4f, %.4f, %.4f] (h, vx, vy)\n', n, X0(1), X0(2), X0(3));
    vx_curr = X0(2);
    [~, idx] = min(abs(vx_range - vx_curr));
    X0_star = X0_stars(:,idx);
    u0_star = u0_stars(:,idx);
    K = K_all{idx};
    fprintf('  K =\n'); disp(K);
    fprintf('  X0_star = [%.4f, %.4f, %.4f], error = [%.4f, %.4f, %.4f]\n', ...
            X0_star(1), X0_star(2), X0_star(3), ...
            X0(1)-X0_star(1), X0(2)-X0_star(2), X0(3)-X0_star(3));
    fprintf('  u0_star = [%.2fº, %.4f, %.4f, %.4f]\n', rad2deg(u0_star(1)), u0_star(2:4));
    fprintf('  K*(X0-X0_star) = [%.4f, %.4f, %.4f, %.4f]\n', (K * (X0 - X0_star))');
    
    u = u0_star + K * (X0 - X0_star);                  % eq 19
    fprintf('  u = [%.2fº, %.4f, %.4f kN/m, %.4f kN/m]\n', rad2deg(u(1)), u(2:4));

    % bound u values to be physically reasonable
    u(1) = max(deg2rad(8), min(deg2rad(30), u(1)));  % th: 8-35 deg
    u(2) = max(deg2rad(-30), min(deg2rad(30), u(2))); % phi: ±30 deg
    u(3) = max(0.5, min(50.0, u(3)));
    u(4) = max(0.5, min(50.0, u(4)));
    
    fprintf('  u clipped = [%.2fº, %.4f, %.4f, %.4f]\n', rad2deg(u(1)), u(2:4));
    % simulate forward one step with adjusted control
    [X1, t_TD, t_LO, com] = slip_return_map(X0, u, params); 
    fprintf('  X1 = [%.4f, %.4f, %.4f]\n', X1);
    
    slip_states(:, n) = X0;
    % add offsets: (OR change all code to allow for pos/t continuity)
    com.p_des(3, :) = com.p_des(3, :) - X0(1);  % convert z to absolute
    pos_log = [pos_log, pos_log(:,end) + com.p_des]; % appending 3 x n_points  
    t_log = [t_log, t_log(end) + com.t_traj];    % appending 1 x n_points
    X0 = X1;
end
save('SLIP3D_data.mat', 'slip_states', 'K', 'X0_star', 'u0_star', 'params', 'pos_log', 't_log');
fprintf('Simulation complete.\n--------------------\n');

%% Functions
% ========================================================================
function [X1, t_TD, t_LO, com] = slip_return_map(X, u, params)
% com: log for ONE sim step N with n_points
%   com.t_traj: time array for entire step (1 x n_points)
%   com.p_des: COM pos traj for entire step (3 x n_points) [x, y, z]
%   com.dp_des: COM vel traj for entire step (3 x n_points) [dx, dy, dz]
%   com.ddp_des: COM accel traj for entire step (3 x n_points) [ddx, ddy, ddz]
% Xs is full slip state [ps, dps]
% X1 is next apex slip state [h;vx;vy]
    
    % currently doesnt allow for pos/t continuity between sim steps N
    h = X(1); vx = X(2); vy = X(3);
    vy = 0; % take out vy change? :) 
    Xs0 = [0; 0; h; vx; vy; 0]; % expand into full SLIP state
    ks1 = u(3);
    ks2 = u(4);
    pf = get_TD_pos(X, u, params);
    tf = params.tf;

    com.t_traj = [];
    com.p_des = [];
    com.dp_des = [];
    com.ddp_des = [];

    % fprintf('θ = %.2f deg\n', rad2deg(u(1))); % debug to check that apex is higher than TD expression lh*cos(th)
    % fprintf('Initial COM height h0 = %.2f\n', Xs0(3));
    % fprintf('Expected TD height (lh*cos(th)) = %.2f\n\n', params.lh*cos(u(1)));
    
    ks = ks1; 
    % Phase 1 - flight: apex to TD
    % change this call so that it doenst start at 0, starts at t0?
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) TDevent(t,Xs,u,params));
    sol = ode45(@(t,Xs) dynamics_SLIP(t,Xs,'flight',params,pf,ks), [0 tf], Xs0, opts);
    if isempty(sol.xe), fprintf('No TD detected during first flight phase'); end
    t0 = sol.x(end);
    start = sol.y(:,end);
    t_TD = sol.xe(end);
    com = update_COM_traj(com, sol, 'flight', params, pf, ks);

    % Phase 2 - stance: TD to max compression (ks1)
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) MCevent(t,Xs,pf));
    sol = ode45(@(t,Xs) dynamics_SLIP(t,Xs,'stance',params,pf,ks), [t0 tf], start, opts);
    if isempty(sol.xe), fprintf('No MC detected during stance phase\n'); end
    t0 = sol.x(end);
    start = sol.y(:,end);
    com = update_COM_traj(com, sol, 'stance', params, pf, ks);
    
    ks = ks2;
    % Phase 3 - stance: max compression to LO (ks2)
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) LOevent(t,Xs,params,pf));
    sol = ode45(@(t,X) dynamics_SLIP(t,X,'stance',params,pf,ks), [t0 tf], start, opts);
    if isempty(sol.xe), fprintf('No LO detected during stance phase\n'); end
    t0 = sol.x(end);
    start = sol.y(:,end);
    t_LO = sol.xe(end);
    com = update_COM_traj(com, sol, 'stance', params, pf, ks);
    
    % Phase 4 - flight: LO to apex
    opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) APevent(t,Xs));
    sol = ode45(@(t,X) dynamics_SLIP(t,X,'flight',params,pf,ks), [t0 tf], start, opts);
    if isempty(sol.xe), fprintf('No apex reached during flight phase\n'); end
    com = update_COM_traj(com, sol, 'flight', params, pf, ks);
    
    % find apex this way to prevent compounding of event detection error???
    [~, idx] = max(sol.y(3,:)); % max z val across all time points
    apex = sol.y(:,idx)';
    % apex = sol.y(:, end);
    X1 = [apex(3); apex(4); apex(5)]; % convert back into simplified apex state [h;vx;vy]
    X1(3) = 0;  % force vy = 0 in output
end

function com = update_COM_traj(com, sol, phase, params, pf, ks)
    t = sol.x; 
    n_points = length(t);

    if strcmp(phase, 'flight')
        ddp = repmat(params.g, 1, n_points); % ballistic, 3xn_points 
    else
        ddp = zeros(3, n_points);
        for i = 1:n_points
            Xs_i = sol.y(:, i);  % state at time point i
            dX = dynamics_SLIP(t(i), Xs_i, 'stance', params, pf, ks);
            ddp(:, i) = dX(4:6);  % only get accel part [ddx; ddy; ddz]
        end
    end
    % update trajectory logs for current sim step N
    % sol.y (6x?) rows = state components, cols = time points
    com.t_traj = [com.t_traj, t];               % concat 1 x n_points
    com.p_des = [com.p_des, sol.y(1:3, :)];     % concat 3 x n_points
    com.dp_des = [com.dp_des, sol.y(4:6, :)];   % concat 3 x n_points
    com.ddp_des = [com.ddp_des, ddp];           % concat 3 x n_points
end

function pf = get_TD_pos(X, u, params)
% eq. 3
% return: foot pos as TD happens
% input: X = [h; vx; vy]. NOT full slip state
    h = X(1);
    theta = u(1);
    phi = u(2);
    % change ps so that it takes in last step ending state
    ps = [0; 0; h];                 % position of mass 3x1 
    phip = [0; -params.yhip; 0];     % position of hip wrt CoM = offset in y-dir 3x1
    params.yhip = -params.yhip;
    th = u(1);
    phi = u(2); 
    lh = params.lh;

    pf = ps + phip + lh * [sin(th)*cos(phi); -sin(th)*sin(phi); -cos(th)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = dynamics_SLIP(~, Xs, phase, params, pf, ks) % time not relevant
% Xs = full slip state
% pf = foot pos at most recent TD (to allow for both flight/stance scenarios)
    M = params.M;
    g = params.g; 
    l0 = params.l0;
    p = Xs(1:3);
    dp = Xs(4:6);
    if strcmp(phase,'flight')   % ballistic dynamics
        dX = [dp; g];
    else                        % stance dynamics: eq. 2
        l = p - pf;
        lhat = l / norm(l); 
        ddp = ks * 1000 * (l0 - norm(l)) * lhat / M + g; % kN/m -> N/m
        dX = [dp; ddp];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVENT DETECTION FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mental note: Xs passed into these must be the full slip state (ps, psd)
function [value, isterminal, direction] = TDevent(~, Xs, u, params) 
% TD event: z pos = l_h * cos(th) = 0 (eq. 3)
    ps = Xs(1:3);
    lh = params.lh;
    th = u(1);
    % fprintf('Xs(3) = %.2f, lh*cos(th) = %.2f\n', Xs(3),params.lh * cos(th));
    value = ps(3) - lh * cos(th); 
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = MCevent(~, Xs, pf)
% MC event: l' * v = 0
    ps = Xs(1:3);
    vs = Xs(4:6);
    l = ps - pf;
    value = l.' * vs;
    isterminal = 1;
    direction = 0;
end

function [value, isterminal, direction] = LOevent(~, Xs, params, pf)
% LO event: ||l|| - l0 = 0, when leg returns to rest length (eq. 4)
    ps = Xs(1:3);
    l = ps - pf;
    value = norm(l) - params.l0;
    isterminal = 1; 
    direction = 1;
end

function [value, isterminal, direction] = APevent(~, Xs) 
% AP event: z vel = 0
    value = Xs(6);
    isterminal = 1;
    direction = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAIT LIBRARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X0_star, u0_star] = find_periodic_gait(X0, params)
% solve LS problem to get (X0*, u0*) optimal state-control pair that
% achieves desired gait timings given a desired forward vel vx
% X0 = [h0; vx; vy0]
% decision variables (z): [h0, vy0, ks, th]

    z0 = [X0(1); X0(3); params.ks0; params.th0]; % initial guess (apex)
    vx = X0(2); % this is kept constant
    fun = @(z) periodic_cost(z, vx, params); 
    options = optimoptions('lsqnonlin','Display','final','MaxFunEvals',2000,'TolX',1e-8);
    % [h0, vy0, ks (kN/m), th]
    lb = [0.21; -1.0; 0.1; deg2rad(15)]; % lh*cos(8º) = 0.79 = min apex height       
    ub = [0.50; 1.0; 30; deg2rad(30)];
    
    % % multi-start: forcing opt to try many initial guesses of th
    % best_sol = [];
    % best_cost = inf;
    % n_starts = 7;  % # of diff starting angles to try
    % th_starts = linspace(deg2rad(15), deg2rad(30), n_starts); % range
    % for i = 1:n_starts
    %     z0 = [X0(1); X0(3); params.ks0; th_starts(i)];
    %     try
    %         sol = lsqnonlin(fun, z0, lb, ub, options);
    %         cost = norm(periodic_cost(sol, vx, params));
    %         if cost < best_cost
    %             best_cost = cost;
    %             best_sol = sol;
    %         end
    %     catch
    %         continue; % skip if opt fails
    %     end
    % end
    % sol = best_sol;
    sol = lsqnonlin(fun, z0, lb, ub, options);
    X0_star = [sol(1); vx; sol(2)];
    u0_star = [sol(4); 0; sol(3); sol(3)];
end

function err = periodic_cost(z, vx, params)
% form symbolic cost function (eq. 13)
    % z = [h0, vy0, ks, th] 4x1 vector of decision vars
    h0  = z(1);
    vy0 = z(2);
    vy0 = 0; % temporarily remove lateral vel from opt
    ks  = z(3);
    th  = z(4);
    X0 = [h0; vx; vy0];     % eq. 14
    u0 = [th; 0; ks; ks];   % eq. 15

    A = diag([1 1 -1]);     % eq. 8, for leg switching (vy)
    [X1, t_TD, t_LO] = slip_return_map(X0, u0, params); % symbolically integrate one step forward
    
    Tdes = get_des_gait_timings(X0); 
    Tcurr = [t_TD; t_LO];   % eq. 12
    
    err = A*X0 - X1;
    % err = [A*X0 - X1; Tdes - Tcurr]; % eq. 13

    % Debug
    % timing_err = abs(Tdes - Tcurr); % Check if timing errors are unreasonably large (infeasible problem)
    % fprintf(    'timing error:\nt_TD: %.4f\nt_LO: %.4f\n', timing_err) % in s

    % make state periodicity (A*X0 = X1) more important
    % weights = [1.0; 1.0; 1.0; 0.01; 0.01];  % [h_err, vx_err, vy_err, t_TD_err, t_LO_err]
    % err = weights .* err(:);

    % Penalty to discourage angles near lower bound (without assuming optimal range)
    % This encourages exploration of higher angles while still allowing lower if truly optimal
    % th = z(4);
    th_lb = deg2rad(15); 
    th_penalty_threshold = deg2rad(2);          % penalty zone: within 2° of lower bound
    if th < th_lb + th_penalty_threshold        % quadratic penalty 
        penalty_weight = 0.05;                  % penalty strength
        th_penalty = penalty_weight * ((th_lb + th_penalty_threshold - th) / th_penalty_threshold)^2;
        err = [err; th_penalty];
    else 
        err = [err; 0];
    end
    % err = err(:); % make sure col vec (can debug w this later)
end

function Tdes = get_des_gait_timings(X0)
% get desired gait timings from human running data to include in LS opt cost
    vx = X0(2);
    c = 2.55*vx^2 - 8.77*vx + 172.9;    % cadence (eq. 9)
    ts = 10^(-0.2) * vx^(-0.82);        % stance time (eq. 11)
    Tstep = 60 / c;                     % period of 1 step 
    Tf = Tstep - ts;                    % flight time = Tstep - ts
    t_TD = Tf / 2;                      % time from apex to TD = half of flight time
    t_LO = Tf / 2 + ts;                 % time to LO = time to TD + stance time
    Tdes = [t_TD; t_LO];
end

function K = compute_deadbeat(X0_star, u0_star, params)
% get gain mat K by using eq. 18 / 19
% where (x0_star, u0_star) is state control pair achieved from LS opt
% Ju du = -Jx dx -> du = K dx 
% so K is -invJu * Jx
    dX = 1e-4; 
    du = 1e-4;
    Jx = zeros(3,3); 
    Ju = zeros(3,4);
    for i = 1:3
        Xp = X0_star; 
        Xp(i)=Xp(i)+dX;
        Xm = X0_star; 
        Xm(i)=Xm(i)-dX;
        Jx(:,i) = (slip_return_map(Xp,u0_star,params) - slip_return_map(Xm,u0_star,params))/(2*dX); % get updated slip state at each perturbation
    end
    for j = 1:4
        up = u0_star; up(j)=up(j)+du;
        um = u0_star; um(j)=um(j)-du;
        Ju(:,j) = (slip_return_map(X0_star,up,params) - slip_return_map(X0_star,um,params))/(2*du);
        % fprintf('Ju column %d: [%.6f, %.6f, %.6f]\n', j, Ju(1,j), Ju(2,j), Ju(3,j));
    end
    
    % du = B * w where w = [dtheta; dphi; dks1]
    B = [1  0  0;
         0  1  0;
         0  0  1;
         0  0 -1];  % dks2 = -dks1

    M = Ju * B;         % Ju * B * w = -Jx * dx, M is a reduced Ju 
    condM = cond(M)
    L = -pinv(M) * Jx;  % w = -inv(M) * Jx * dx, so L = -pinv(M) * Jx

    % eps_reg = 1e-4;  % Start with this for cond(M) ~ 4000
    % L = -pinv(M + eps_reg * eye(3)) * Jx;

    K = B * L; % K = B * -pinv(M) * Jx, ks1 row and ks2 row will be negatives

    % K = -pinv(Ju)*Jx; % original method 
end