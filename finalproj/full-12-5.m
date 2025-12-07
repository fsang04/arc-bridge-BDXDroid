%% 12/4-12/5-12/6 gait library adding LCM
% copy of gaitlib-11-27.m
% using draft from slip3D_WBC_integrated.m

clear; clc;
run ../setup.m 

%% LCM Setup
global lcm_state_topic lcm_cmd_topic
lcm_state_topic = "bdx_droid_state";
lcm_cmd_topic   = "bdx_droid_control";

lc = lcm.lcm.LCM.getSingleton();
getenv("LCM_DEFAULT_URL")

aggregator = lcm.lcm.MessageAggregator();
aggregator.setMaxMessages(1);
lc.subscribe(lcm_state_topic, aggregator);

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
params.mu = 1; % friction
params.dt = dt;

% sim parameters
N = 10; % number of steps
params.X0 = [0.21; 0.6; 0]; % starting X for sim loop

% % robot physical parameters (for dynamics computation)
% params.l = [0.2; 0.0; 0.4; 0.4; 0.03];  % [l_hip_roll; l_hip_pitch; l_thigh; l_shin; l_foot]
% params.M_masses = [1.685713; 0.158036; 0.025867; 0.026811];  % [M_trunk; M_thigh; M_shin; M_foot]
% params.I_inertias = [0.01; 0.005; 0.005; 7.0e-05];  % [I_trunk; I_thigh; I_shin; I_foot]
% params.params = [params.g(3); params.l; params.M_masses; params.I_inertias];

%% Build gait library 
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

%% Simulation test (within Matlab no LCM)
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
    % logging is currently non-cumulative across steps
    % add offsets: (OR change all code to allow for pos/t continuity)
    com.p_des(3, :) = com.p_des(3, :) - X0(1);          % get rid of z-offset from last
    pos_log = [pos_log, pos_log(:,end) + com.p_des];    % appending 3 x n_points  
    t_log = [t_log, t_log(end) + com.t_traj];           % appending 1 x n_points
    X0 = X1;
end
save('SLIP3D_data.mat', 'slip_states', 'K', 'X0_star', 'u0_star', 'params', 'pos_log', 't_log');
fprintf('Simulation complete.\n--------------------\n');

%% WBC Parameters
Kp_com = diag([1000, 1000, 1500]);  % pos gains [x, y, z]
Kd_com = diag([100, 100, 150]);     % vel gains 
% foot controller gains (4x4 for both feet: [Rx; Rz; Lx; Lz])
params.Kp_p = diag(repmat([500, 800], [1, 2]));     % pos gains [Rx, Rz, Lx, Lz]
params.Kd_p = diag(repmat([50, 60], [1, 2]));       % vel gains [Rx, Rz, Lx, Lz]
% params.Kp_w;
% params.Kd_w;
params.w_tau = 1;
params.w_qdd = 1e-5;
params.w_Fs = 1e-4;

% task weighting: 4 (foot) + 3 (com) + 3 (ang mom) + 10 (pose)
w_foot = ones(1, 4);                  % [1, 1, 1, 1]
w_com = repmat(25, [1, 3]);           % [25, 25, 25]
w_angmom = [20, 4, 20];               % [20, 4, 20]
% w_torso = [17.5, 70, 14];             % [17.5, 70, 14]
w_hip = repmat(0.1, [1, 6]);          % [0.1, 0.1]
w_knee = repmat(0.5, [1, 2]);         % [0.5, 0.5]
w_ankle = repmat(0.1, [1, 2]);        % [0.1, 0.1] 
w_task_vec = [w_foot, w_com, w_angmom, w_hip, w_knee, w_ankle];
params.w_task = diag(w_task_vec);     % 20x20 diagonal weighting matrix

rs = struct();
rs.Xslip = params.X0;               % will store current SLIP apex state [h; vx; vy]
rs.uslip = [deg2rad(25); 0; 1500; 1500]; % will store SLIP control [th; phi; ks1; ks2]
rs.stance_foot_idx = 1;             % 1 = left foot stance, 2 = right foot stance
t = 0;                           % keep track of curr time
rs.step_count = 0;                  % step counter

import casadi.*
PTSC_func = setup_PTSC(params);

% %% LCM Main Control Loop
fprintf('Starting TSC control loop...\n');
while true
    msg = aggregator.getNextMessage(0);
    if isempty(msg)
        continue;
    end
    lc_state = eval("lcm_msgs."+lcm_state_topic+"_t(msg.data)");
    lc_cmd   = eval("lcm_msgs."+lcm_cmd_topic+"_t()");
    
    % ------- 1. SLIP template update -------
    rs = updateRobotSLIPState(t, lc_state, rs, params); 

    % retrieve closest K corresponding to vx from "library":
    vx_curr = rs.Xslip(2);
    [~, idx] = min(abs(vx_range - vx_curr));
    X0_star = X0_stars(:, idx);
    u0_star = u0_stars(:, idx);
    K = K_all{idx};
    fprintf('Optimal library gait found.\n');
    fprintf('  K =\n'); disp(K);
    fprintf('  rs.Xslip = [%.2fº, %.4f, %.4f, %.4f]\n', rs.Xslip);
    fprintf('  X0_star = [%.4f, %.4f, %.4f], error = [%.4f, %.4f, %.4f]\n', ...
            X0_star(1), X0_star(2), X0_star(3), ...
            X0(1)-X0_star(1), X0(2)-X0_star(2), X0(3)-X0_star(3));
    fprintf('  u0_star = [%.2fº, %.4f, %.4f, %.4f]\n', rad2deg(u0_star(1)), u0_star(2:4));
    fprintf('  K*(X0-X0_star) = [%.4f, %.4f, %.4f, %.4f]\n', (K * (X0 - X0_star))');
    u = u0_star + K * (rs.Xslip - X0_star);                  % eq 19
    fprintf('  u = [%.2fº, %.4f, %.4f kN/m, %.4f kN/m]\n', rad2deg(u(1)), u(2:4));
    rs.uslip = u;   % store to be sent over later

    % bound u values to be physically reasonable
    u(1) = max(deg2rad(8), min(deg2rad(30), u(1)));  % th: 8-30 deg
    u(2) = max(deg2rad(-30), min(deg2rad(30), u(2))); % phi: ±30 deg
    u(3) = max(0.5, min(50.0, u(3)));
    u(4) = max(0.5, min(50.0, u(4)));

    % ------- 2. update desired CoM traj / foot traj -------
    [X1, t_TD, t_LO, com, pf_TD] = slip_return_map(rs.Xslip, rs.uslip, params);
    rs.t_TD = t_TD; % can integrate rs into the function later 
    rs.t_LO = t_LO;
    rs.pf_TD = pf_TD;
    rs.T_step = com.t_traj(end);  % total step period (bc time is non-cumulative)
    
    % store pf_LO at LO time (used to be in updateRobotSLIPState, moved bc it only happens once per loop)
    t_step = mod(t - rs.t_step_start, rs.T_step);
    dt = 1/500;
    if abs(t_step - rs.t_LO) < dt % check if within one timestep of t_LO -> prone to error?
        rs.pf_trans(:, 1) = rs.pf(1:2);  % right foot lifts off. store right foot pos [x; z]
    end

    % ------- 3. compute prioritized tasks ------
    rsDes = updateRobotStateDes(t, rs);
    [Atask, btask] = compute_prioritized_tasks(rs, rsDes, params);
    
    % ------- 4. PTSC -------
    [tau, qdd, Fs] = my_PTSC(PTSC_func, rs, Atask, btask);
 
    % ------- 5. publish to LCM -------
    lc_cmd.qj_tau = tau;  % (10 x 1) actuated joint torques
    lc.publish(lcm_cmd_topic, lc_cmd);

    rate_ctrl.waitfor();
    t = t + dt;
end

% ========================================================================
% %% LCM functions
% ========================================================================
function rs = updateRobotSLIPState(t, state, rs, p)
    pos = state.position;
    vel = state.velocity;
    rpy = state.rpy;
    omega = state.omega;
    quat = state.quaternion;
    qj = state.qj_pos;
    dqj = state.qj_vel;
    pf = state.p_gc;            % foot positions wrt world (task space) [right x; right z; left x; left z]
    Jf = state.J_gc;
    dJfdq = state.dJdq_gc;
    ptor = state.position;      % is this correct? 
    Jtor = state.J_tor;
    dJtordq = state.dJdq_tor;

    % NOTE: either use LCM or matlab FK functions?
    % [ptor, Jtor, dJtordq] = fcn_torso_p_J_dJdq(rs.q, rs.dq, p.l);
    % [H, bias] = fcn_droid_Mass_bias(rs.q, rs.dq, p.params);
    H = state.inertia_mat;      % 16x16 -> convert to 13x13 ?
    bias = state.bias_force;    % 16x1 -> convert to 13x1 ?
    % % conversion example code (16x16 -> 13x13)
    % % LCM: 6 floating (x,y,z,roll,pitch,yaw) + 10 actuated
    % % We use: 3 floating (x,z,pitch) + 10 actuated
    % H_lcm = reshape(lc_state.inertia_mat, 16, 16);
    % bias_lcm = lc_state.bias_force(:);
    % 
    % % Map LCM DoF to our DoF: [x(1), z(3), pitch(5), q4(7), ..., q13(16)]
    % lcm_to_our_dof = [1, 3, 5, 7:16];  % x, z, pitch, then all actuated joints
    % H = H_lcm(lcm_to_our_dof, lcm_to_our_dof);
    % bias = bias_lcm(lcm_to_our_dof);
    
    rs.qb = [pos(1);pos(3);rpy(2)];     % base [x,z,th]
    rs.dqb = [vel(1);vel(3);omega(2)];
    rs.acc = state.acceleration([1,3]);
    rs.qj = qj;                         % actuated joints
    rs.dqj = dqj;
    rs.q = [rs.qb; rs.qj];              
    rs.dq = [rs.dqb; rs.dqj];
    rs.Xslip = [pos(3); vel(1); vel(2)]; % Xslip = [h; vx; vy]

    rs.ptor = ptor;
    rs.Jtor = Jtor;         % just edited bdx_droid_bridge.py
    rs.vtor = Jtor * rs.dq; % torso vel wrt slip CoM
    rs.dJtordq = dJtordq;   % also edited bdx bridge
    
    % 2d control for now
    pf_reshape = [pf(1); pf(3); pf(4); pf(6)]; % 6x1 -> 4x1
    Jf_reshape = [Jf(1,:); Jf(3,:); Jf(4,:); Jf(6,:)]; % 6x13 -> 4x13
    dJfdq_reshape = [dJfdq(1); dJfdq(3); dJfdq(4); dJfdq(6)];
    rs.pf = pf_reshape;              
    rs.Jf = Jf_reshape;
    rs.vf = rs.Jf * rs.dq;     % foot vel wrt slip CoM (task space)
    rs.dJfdq = dJfdq_reshape;
    rs.H = H;
    rs.bias = bias;
    
    % init at start of sim
    if t < p.dt
        rs.t_step_start = t;
        rs.stance_foot_idx = 1; % start with right foot in stance          
        rs.pf_trans = reshape(rs.pf, [2,2]);    % [right; left] columns
        fprintf('rs.pf_trans = \n'); disp(rs.pf_trans);
    end
    % update stance foot and detect LO
    if isfield(rs, 't_TD') && isfield(rs, 't_LO') && isfield(rs, 'T_step')
        % check current time in current step cycle
        t_step = mod(t - rs.t_step_start, rs.T_step); % [0, T_step)
        
        % right leg based tracking: if in stance phase (t_TD to t_LO), right leg is stance
        if (t_step >= rs.t_TD) && (t_step < rs.t_LO)    % stance: t = [t_TD, t_LO] 
            rs.stance_foot_idx = 1;                     % right leg stance
        else
            rs.stance_foot_idx = 2;                     % left leg
        end
    end
end

function rsDes = updateRobotStateDes(t, rs) 
    t_LO = rs.t_LO;
    p_TD = rs.pf_TD; % stored from main loop
    t_step = mod(t - rs.t_step_start, rs.T_step);  % time within current step
    t_swing = t - t_LO;
   
    pfd = zeros(4, 1); % [Rx; Rz; Lx; Lz]
    dpfd = zeros(4, 1);
    ddpfd = zeros(4, 1);

    % NOTE: migt need to fix for x;y;z instead of x;z?
    % determine if each leg is in stance
    for i_leg = 1:2  % 1 = R, 2 = L
        idx = (i_leg-1) * 2 + (1:2);
        if rs.stance_foot_idx == 1   % right stance
            is_swing = (i_leg == 2); % in swing if idx = left leg
        else 
            is_swing = (i_leg == 1);
        end
        if is_swing % leg is currently in swing
            p_LO = rs.pf_trans(:, i_leg);  % get LO pos [x; z]
            p_TD = rs.pf_TD;  % get TD pos
    
            % compute how long leg has been in swing 
            if t_step < rs.t_TD % flight phase 1 (AP -> TD)
                t_swing = t_step;               % time since AP
                T_flight = rs.t_TD;             % duration of first flight
            else                % flight phase 2 (LO -> AP)
                t_swing = t_step - rs.t_LO;     % time since LO
                T_flight = rs.T_step - rs.t_LO; % duration of second flight
            end
            % normalize swing progress [0, 1]
            s_sw = t_swing / T_flight;
            s_sw = max(0, min(s_sw, 1.0));  % clamp
    
            [pf_des, dpf_des, ddpf_des] = cubic_spline(p_TD, p_LO, s_sw, T_flight); % get desired foot values at current time idx
            pfd(idx) = pf_des;  % populate corresponding part of pfd with a pos [x; z]
            dpfd(idx) = dpf_des; 
            ddpfd(idx) = ddpf_des;
        else % leg is currently in stance, keep foot fixed
            pfd(idx) = rs.pf(idx);  % current pos
            dpfd(idx) = zeros(2, 1);
            ddpfd(idx) = zeros(2, 1);
        end
    end

    rsDes.pf = pfd;
    rsDes.dpf = dpfd;
    rsDes.ddpf = ddpfd;

    % compute w, w_des, wd_des, 
end

function [p, v, a] = cubic_spline(p0, pf, s, T)
% s: normalized time [0,1]
% T: flight duration
% outputs: at a single idx
%   p: desired swing foot pos (2 x 1) [x; z]
%   v: desired swing foot vel (2 x 1) [dx; dz]
%   a: desired swing foot accel (2 x 1) [ddx; ddz]
    % boundary conditions
    v0 = [0; 0];  % init velocity (can be from actual foot vel?)
    vf = [0; 0];  % final vel (zero as defined in paper)
    
    % cubic polynomial coefficients (T = 1.0)
    a0 = p0;
    a1 = v0 * T;
    a2 = 3*(pf - p0) - (2*v0 + vf)*T;
    a3 = -2*(pf - p0) + (v0 + vf)*T;
    
    % evaluate polynomial
    p = a0 + a1*s + a2*s^2 + a3*s^3;
    v = (1/T) * (a1 + 2*a2*s + 3*a3*s^2); % dp/ds -> dp/dt
    a = (1/T^2) * (2*a2 + 6*a3*s); % ddp/ds -> dp/dt

    % add ground clearance (parabolic lift)
    clearance = 0.05; % could adjust
    z_clear = 4 * clearance * s * (1 - s); % clearance pos
    dz_clear = 4 * clearance * (1 - 2*s) / T; % clearance vel
    ddz_clear = -8 * clearance / T^2; % clearance accel
    
    % Apply clearance to z-component only
    p(2) = p(2) + z_clear;
    v(2) = v(2) + dz_clear;
    a(2) = a(2) + ddz_clear;
end

% ========================================================================
% %% TSC functions
% ========================================================================
function e_theta = orientationError(wd, w)
    % e_theta: 3x1 angle-axis representation of error between a desired and actual orientation
end

function [ddp_c, dw_c] = foot_controller(rs, rsDes, params)    % eq. 22/23
% outputs
%   ddp_c: commanded linear foot accel 4x1 [Rx; Rz; Lx; Lz]
%   wd_c: commanded angular foot accel 
% if in stance: ddp_c = 0 + Kd*(0 - dpf_act) + Kp*(pf_act - pf_act) = -Kd*dpf_act
% paper directly sets to 0, current code does not
    % e_theta = orientationError(w_des, w);
    % dw_c = wd_des + Kd_w * (w_des - w) + Kp_w * e_theta;
    dw_c = zeros(2, 1); % for now 
    ddp_c = rsDes.ddpf + params.Kd_p * (rsDes.dpf - rs.vf) + params.Kp_p * (rsDes.pf - rs.pf);
end

function [dl_Gc, dk_Gc] = momentum_controller(com, rs, params, Kp_l, Kd_l, Kd_k)    % eq. 24/25
%   ps, dps, ddps: SLIP CoM pos, vel, accel (desired)
%   pG, dpG: robot CoM pos, vel (actual)
%   Kp_l, Kd_l: gains
%   Kd_k: gains
% outputs: 
%   ld_Gc: commanded rate of change in total system linear momentum (from
%   PD control of humanoid CoM to 3D slip)
%   kd_Gc: commanded rate of change in centroidal angular momentum

    ps = com.p_des(:, current_idx); 
    dps = com.dp_des(:, current_idx); 
    ddps = com.ddp_des(:, current_idx); 
    pG = rs.ptor(1:3);  
    dpG = rs.vtor(1:3);
    m = params.m;
    
    % [A, Ad_qd] = CMM(model, rs.q, rs.dq); % A = 6xN CMM
    h = A * rs.dq; % 6 x 1 
    kG = h(4:6);  % centroidal angular momentum (3 x 1)
    
    % For 3D case, k_G would be (3 x 1)
    % For 2D case, k_G is scalar (pitch component only)

    dl_Gc = m*(ddps + Kd_l*(dps - dpG) + Kp_l*(ps - pG));
    dk_Gc = -Kd_k * kG; 
end

function ddq_c = pose_controller(rs, params)
% for actuated joints only
% dqj_des, ddqj_des = 0 (from paper)
% ddq_c is 10x1
    qj_des = params.q_nominal;  % (10 x 1) fixed desired poses -> get qj_des from slip kinematics? idk
    dqj_des = zeros(N_actuated, 1);  % zero desired vel
    ddqj_des = zeros(N_actuated, 1);  % zero desired accel
    
    Kp = params.KP_pose;  % (10 x 1) or scalar
    Kd = params.KD_pose;  % (10 x 1) or scalar

    ddq_c = ddqj_des + Kd*(dqj_des - rs.dqj) + Kp*(qj_des - rs.qj);  % eq. 26, revolute joints
end

function [Atask, btask] = compute_prioritized_tasks(rs, rsDes, params)
% need:
%   ddp_c, dw_c: commanded foot controller dynamics
%   dl_Gc, dk_Gc: commanded momentum controller dynamics
%   ddq_c: commanded pose controller dynamics
    
    % FOOT TASK:
    % J_foot * qdd + dJfdq = ddp_c -> J_foot * qdd = ddp_c - dJfdq
    [ddp_c, dw_c] = foot_controller(rs, rsDes, params); 
    A1 = rs.Jf; % foot jacobian 6x13 -> change to 4x13 [Rx; Rz; Lx; Lz] for now 
    b1 = ddp_c - rs.dJfdq;  % 4x1

    % will add null space projection when adding other tasks
    % A2 = J_mom * N1, where N1 = null space of A1
    % A3 = J_pose * N2, where N2 = null space of [A1; A2]
    % MOMENTUM TASK: (TODO)
    A2 = zeros(6, 13); 
    b2 = zeros(6, 1);

    % POSE TASK: (TODO)
    A3 = zeros(10, 13); 
    b3 = zeros(10, 1);

    Atask = [A1; A2; A3]; % 20x13
    btask = [b1; b2; b3]; % 20x1
end

function [tau, qdd, Fs] = my_PTSC(PTSC_func, rs, Atask, btask)
% call compiled casadi func
import casadi.* % to convert to DM
    contact = double([rs.stance_foot_idx == 1; rs.stance_foot_idx == 2]);  % 2x1
    [tau, qdd, Fs] = PTSC_func(rs.q, rs.dq, Atask, btask, rs.H, rs.bias, rs.Jf, contact); 
    tau = full(tau);
    qdd = full(qdd);
    Fs = full(Fs);
end

function PTSC_func = setup_PTSC(params) 
% compile casadi symbolic func 
% params passed at runtime (values that are updated every time PTSC is called):
%   q, qd: nx1
%   tasks: Atask task jacobian, btask task error
%   H: nxn mass matrix
%   bias: nX1 bias Cqdot + G
%   J_c: combined contact Jacobian (Jf from bdx lcm)
%   contact: contact flags [right_foot, left_foot]
    
    N = 13;     % num DoF
    Ntau = 10;  % num actuated DoF
    Ntask = 20; % 4 (foot) + 6 (mom) + 10 (pose)
    nc = 4;     % num contact constraints (2 per foot?)
    mu = params.mu;
    
    import casadi.*
    opti = casadi.Opti('conic');

    % decision vars: [tau, F_c, ddq]
    tau = opti.variable(Ntau, 1);   % joint torques 10x1
    qdd = opti.variable(N, 1);      % joint accels 13x1
    Fs = opti.variable(4, 1);       % 4x1 for 2D contact (x,y)

    % parameters (passed at runtime):
    q = opti.parameter(N, 1);
    dq = opti.parameter(N, 1);
    Atask = opti.parameter(Ntask, N);  % total task Jacobian (with priorities)
    btask = opti.parameter(Ntask, 1);  % total task error (with priorities)
    H = opti.parameter(N, N);
    bias = opti.parameter(N, 1);
    J_c = opti.parameter(nc, N); 
    contact = opti.parameter(2, 1);   

    % contact selection matrix
    Sa = [zeros(Ntau,3), ones(Ntau,Ntau)]; % actuation selection matrix [0_10x6 1_10x10] 3 or 6?

    % objective: minimize task tracking error + regularization
    e_task = Atask * qdd + btask;  % task error
    obj = MX(0);
    obj = obj + e_task' * params.w_task * e_task;  % ||Atask*qdd + btask||^2
    obj = obj + tau' * params.w_tau * tau;    % weighted torques
    obj = obj + qdd' * params.w_qdd * qdd;    % weighted accelerations
    obj = obj + Fs' * params.w_Fs * Fs; % weighted contact 
    opti.minimize(obj);
    
    % constraints: 
    opti.subject_to(H * qdd + bias == Sa' * tau + J_c' * Fs); % dynamics constraint
    F_R = Fs(1:2);
    F_L = Fs(3:4);
    opti.subject_to(-mu * F_R(2) <= F_R(1) <= mu * F_R(2)); % friction cone constraints: -mu*Fz <= Fx <= mu*Fz
    opti.subject_to(-mu * F_L(2) <= F_L(1) <= mu * F_L(2));
    opti.subject_to(F_R(2) <= contact(1) * 500 + 0.5); % max force when in contact 
    opti.subject_to(F_L(2) <= contact(2) * 500 + 0.5);
   
    opti.solver('osqp');
    input = {q, dq, Atask, btask, H, bias, J_c, contact};
    output = {tau, qdd, Fs};

    PTSC_func = opti.to_function('prioritized_TSC', input, output);
end

% ========================================================================
%% Functions
% ========================================================================
function [X1, t_TD, t_LO, com, pf_TD] = slip_return_map(X, u, params)
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
    pf_TD = pf([1, 3]); % 2d ? 
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
    % update trajectory logs for CURRENT sim step N (non-cumulative)
    % sol.y (6x?) rows = state components, cols = time points
    com.t_traj = [com.t_traj, t];               % concat 1 x n_points
    com.p_des = [com.p_des, sol.y(1:3, :)];     % concat 3 x n_points
    com.dp_des = [com.dp_des, sol.y(4:6, :)];   % concat 3 x n_points
    com.ddp_des = [com.ddp_des, ddp];           % concat 3 x n_points
end

function pf = get_TD_pos(X, u, params)
% eq. 3
% return: foot pos of next TD [x; y; z] 
% input: X = [h; vx; vy]. NOT full slip state
    h = X(1);
    theta = u(1);
    phi = u(2);
    % change ps so that it takes in last step ending state
    ps = [0; 0; h];                 % position of mass 3x1 
    phip = [0; -params.yhip; 0];     % position of hip wrt CoM = offset in y-dir 3x1
    
    % % hip offset sign depends on which leg is in stance?
    % if stance_foot_idx == 1  % left leg stance
    %     phip = [0; params.yhip; 0];
    % else  % right leg stance
    %     phip = [0; -params.yhip; 0];
    % end
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

    % debug
    % timing_err = abs(Tdes - Tcurr); % check if timing errors are too large
    % fprintf(    'timing error:\nt_TD: %.4f\nt_LO: %.4f\n', timing_err) % in s

    % make state periodicity (A*X0 = X1) more important
    % weights = [1.0; 1.0; 1.0; 0.01; 0.01];  % [h_err, vx_err, vy_err, t_TD_err, t_LO_err]
    % err = weights .* err(:);

    % penalty to discourage angles near lower bound (without assuming optimal range)
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