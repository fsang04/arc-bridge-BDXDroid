%% slip3D with LCM simulation loop, no TSC
%% minimum LCM control needed to interface with MuJoCo
clear; clc;
run ../setup.m
 
%% hi
FK_FD_droid()

%% LCM Setup
global lcm_state_topic lcm_cmd_topic
lcm_state_topic = "bdx_droid_state";
lcm_cmd_topic   = "bdx_droid_control";

lc = lcm.lcm.LCM.getSingleton();
getenv("LCM_DEFAULT_URL")

aggregator = lcm.lcm.MessageAggregator();
aggregator.setMaxMessages(1);
lc.subscribe(lcm_state_topic, aggregator);

%% Control Parameters
control_freq = 200;  % Control frequency in Hz (200 in paper)
rate_ctrl = rateControl(control_freq);
dt = 1 / control_freq;

%% SLIP Parameters
params.M = 5;               % effective point mass
params.g = [0; 0; -9.81];
params.l0 = 0.24;           % rest spring leg length, at TD l0 = lh
params.lh = 0.24;            % humanoid virtual leg length used to map to SLIP leg
params.yhip = 0.0;        % zero for testing. hip offset in y-dir (left hip at y=0.035, right hip at y=-0.035)
params.th0 = deg2rad(24);   % init TD angle guess
params.ks0 = 6000;          % init stiffness guess
params.tf = 5.0;            % single step time interval

% robot physical parameters (for dynamics computation)
params.l = [0.2; 0.0; 0.4; 0.4; 0.03];  % [l_hip_roll; l_hip_pitch; l_thigh; l_shin; l_foot]
params.M_masses = [10.0; 1.0; 0.5; 0.05];  % [M_trunk; M_thigh; M_shin; M_foot]
params.I_inertias = [0.01; 0.005; 0.005; 7.0e-05];  % [I_trunk; I_thigh; I_shin; I_foot]
params.params = [params.g(3); params.l; params.M_masses; params.I_inertias];

%% Load SLIP gait library computed offline
gait_library_file = 'SLIP3D_gait_library.mat';
if exist(gait_library_file, 'file')
    fprintf('Loading SLIP gait library from %s...\n', gait_library_file);
    load(gait_library_file, 'vx_range', 'X0_stars', 'u0_stars', 'K_all');
    fprintf('Gait library loaded: %d gaits\n', length(vx_range));
else
    fprintf('Warning: Gait library not found. Will use default SLIP parameters.\n');
    vx_range = [];
    X0_stars = [];
    u0_stars = [];
    K_all = {};
end

rs = struct();
rs.Xslip = [];      % will store current SLIP apex state [h; vx; vy]
rs.uslip = [];    % will store SLIP control [theta; phi; ks1; ks2]
rs.stance_foot_idx = 1;  % 1 = left foot stance, 2 = right foot stance
rs.t = 0;                
rs.step_count = 0;       

% init SLIP state (will be updated from robot state)
X_slip = [2.0; 3.5; 0.0];  % [h; vx; vy] - initial apex state
u_slip = [deg2rad(24); 0; 6000; 6000];  % [theta; phi; ks1; ks2]

t = 0;

% control gains (minimal PD control)
Kp_com = 10000 * eye(3);  % pos gain for CoM
Kd_com = 100 * eye(3);   % vel gain for CoM

%% Main Loop
while true
    msg = aggregator.getNextMessage(0);
    if isempty(msg)
        rate_ctrl.waitfor();
        rs.t = rs.t + dt;
        continue;
    end

    lc_state = eval("lcm_msgs."+lcm_state_topic+"_t(msg.data)");
    lc_cmd   = eval("lcm_msgs."+lcm_cmd_topic+"_t()");
    
    % ------- 1. parse robot state from LCM -------
    % extract basic state information
    rs = updateRobotSLIPState(t, lc_state, rs);
    
    % current CoM state
    x_com = rs.qb;  % [x, y, z] (3x1)
    
    % get CoM jac from torso jac
    % Jtor is 6x13: [linear_vel (3x13); angular_vel (3x13)]
    J_com = rs.Jtor(1:3, :);  % 3x13 - extract linear velocity rows
    
    % compute CoM velocity in task space (3x1)
    xd_com = J_com * rs.dq;  % (3x13) @ (13x1) = (3x1) -> [vx, vy, vz]
    
    % ------- 2. update SLIP state and compute control -------
    % = [h; vx; vy]
    Xslip = rs.Xslip;
    
    % get SLIP control from gait library
    if ~isempty(vx_range) && length(vx_range) > 0
        fprintf('Retrieving from library...\n');
        fprintf('Step: Xcom = [%.4f, %.4f, %.4f] (x, y, z)\n', x_com(1), x_com(2), x_com(3));
        fprintf('Step: Xslip = [%.4f, %.4f, %.4f] (h, vx, vy)\n', Xslip(1), Xslip(2), Xslip(3));
        vx_curr = Xslip(2);
        [~, idx] = min(abs(vx_range - vx_curr));
        if idx <= length(vx_range)
            X0_star = X0_stars(:, idx);
            u0_star = u0_stars(:, idx);
            K = K_all{idx};
            
            u = u0_star + K * (Xslip - X0_star);  % deadbeat control (eq. 19)
        else
            u = [deg2rad(24); 0; 6000; 6000];  % default
        end
    else
        u = [deg2rad(24); 0; 6000; 6000];  % default if no library
    end
    
    % bound u values to be physically reasonable
    u(1) = max(deg2rad(8), min(deg2rad(35), u(1)));  % th: 8-35 deg
    u(2) = max(deg2rad(-30), min(deg2rad(30), u(2))); % phi: ±30 deg
    u(3) = max(1000, min(50000, u(3))); 
    u(4) = max(1000, min(50000, u(4)));
    
    % ------- 3. minimal PD control -------
    % convert SLIP control to desired CoM position
    % desired CoM: maintain height from SLIP, track forward velocity
    x_com_des = [rs.qb(1) + rs.dqb(1) * dt;
             rs.qb(2);   % maintain current y position (passive), this should stay near 0
             Xslip(1)];   % z from SLIP
    
    xd_com_des = [rs.dqb(1);  % vx from SLIP
              0;  % currently constrain y vel to 0
              rs.dqb(1)];         % vz = 0 (no vertical motion at desired)
    
    e_pos = x_com_des - x_com;
    e_vel = xd_com_des - xd_com;
    F_com = Kp_com * e_pos + Kd_com * e_vel;  % task space force

    tau_full = J_com' * F_com;  % 13x1 (includes floating base + actuated)
    
    % actuated joint torques last 10 DoF: q4-q13
    tau_actuated = tau_full(4:13);  % 10x1
    
    % ------- 4. publish LCM control -------
    lc_cmd.timestamp = java.lang.System.nanoTime();
    lc_cmd.qj_tau = tau_actuated(:);  % 10x1 joint torques
    lc_cmd.qj_pos = zeros(10, 1);     % not using position control
    lc_cmd.qj_vel = zeros(10, 1);     % not using velocity control
    lc_cmd.kp = zeros(10, 1);         % no PD gains
    lc_cmd.kd = zeros(10, 1);
    lc_cmd.reset_se = false;
    lc_cmd.se_pos = zeros(3, 1);
    lc_cmd.se_vel = zeros(3, 1);
    lc_cmd.EA = zeros(3, 1);
    
    lc.publish(char(lcm_cmd_topic), lc_cmd);
    
    rate_ctrl.waitfor();
    t = t + dt;
end

function rs = updateRobotSLIPState(t, state, rs)
    pos = state.position; % [x, y, z]
    vel = state.velocity; % [vx, vy, vz]
    rpy = state.rpy; % [roll, pitch, yaw]
    omega = state.omega; %  [wx, wy, wz]
    quat = state.quaternion;
    qj = state.qj_pos;  % 13x1 joint positions
    dqj = state.qj_vel; % 13x1 joint velocities
    pf = state.p_gc;
    Jf = state.J_gc;
    dJfdq = state.dJdq_gc;
    ptor = state.position;      % is this correct? 
    Jtor = state.J_tor; % 
    dJtordq = state.dJdq_tor; %
    H = state.inertia_mat;      % 13x13 
    bias = state.bias_force; 

    % rs.qb = [pos(1);pos(3);rpy(2)];     % [x,z,th]
    % rs.dqb = [vel(1);vel(3);omega(2)];
    rs.qb = [pos(1);pos(2);pos(3)];
    rs.dqb = [vel(1);vel(2);vel(3)];
    rs.acc = state.acceleration([1,3]);
    rs.qj = qj;
    rs.dqj = dqj;
    rs.q = [rs.qb; rs.qj];
    rs.dq = [rs.dqb; rs.dqj];
    rs.Xslip = [pos(3); vel(1); vel(2)]; % Xslip = [h; vx; vy]

    rs.ptor = ptor;
    rs.Jtor = Jtor;          % just edited bdx_droid_bridge.py
    rs.vtor = Jtor * rs.dq;
    rs.dJtordq = dJtordq;    % also edited bdx bridge
    rs.pf = pf;
    rs.Jf = Jf;
    rs.vf = Jf * rs.dq;
    rs.dJfdq = dJfdq;
    rs.H = H;
    rs.bias = bias;
end