%% draft v2
%% main loop

clear; clc;
run ../setup.m 
%%
FK_FD_droid()
%%

global lcm_state_topic lcm_cmd_topic
lcm_state_topic = "bdx_droid_state";
lcm_cmd_topic   = "bdx_droid_control";

lc = lcm.lcm.LCM.getSingleton();
getenv("LCM_DEFAULT_URL")

aggregator = lcm.lcm.MessageAggregator();
aggregator.setMaxMessages(1);
lc.subscribe(lcm_state_topic, aggregator);

control_freq = 500; % control frequency in Hz
rate_ctrl = rateControl(control_freq);
dt = 1 / control_freq;
steps = 1000; % Tot steps for the simulation

p.N     = 50;          % collocation intervals
p.tf    = 0.6;         % total duration of jump [s], tune this
p.g     = 9.81;         
p.z_min = 0.4;         % minimum CoM height

% 
% p.l: [l_hip_roll; l_hip_pitch; l_thigh; l_shin; l_foot]
p.l = [0.2; 0.0; 0.4; 0.4; 0.03];
% p.M: [M_trunk; M_thigh; M_shin; M_foot]
p.M = [10.0; 1.0; 0.5; 0.05];
% p.I: [I_trunk; I_thigh; I_shin; I_foot]
p.I = [0.01; 0.005; 0.005; 7.0e-05];

% sum M for SLIP cart model
p.M_eff = sum(p.M); % should be M(1)?

% Parameters vector for mass/bias function
p.params = [p.g; p.l; p.M; p.I];

% SLIP leg parameters (spring-damper model)
p.k  = 5000;   % stiffness 
p.c  = 200;    % damping 
p.l0 = 0.6;
p.Umax  = [3000; 4000];  % [Fx_max; Fz_max], tune?

% stance vs flight: both feet together for synchronized jump
stance_frac = 0.4;  % 40% of motion is stance phase
Ns          = floor(stance_frac * p.N);
% contact_mask is 1x(N+1): 1 = both feet in contact, 0 = both feet in flight
p.contact_mask = [ones(1,Ns+1), zeros(1, p.N-Ns)];

% objective weights
p.w_dist = 200;    % weight for jump distance 
p.w_vf   = 5;      % weight for final velocity
p.use_slip = true; % tells updateRobotState not to override ifStance

% WBC-related params you already have:
p.mu     = 0.5;
p.alpha  = 0;
p.Q_tor  = [1 1 1];
p.Q_f    = [1 1 1 1];

controller_WBC = ctrl_WBC(p);
rs = initRobotState(p);

jump_traj = slip_traj_opt(p, rs);

lc = lcm.lcm.LCM.getSingleton();
aggregator = lcm.lcm.MessageAggregator();
aggregator.setMaxMessages(1);
lc.subscribe(char(lcm_state_topic), aggregator);

control_freq = 500;
rate_ctrl    = rateControl(control_freq);
dt           = 1 / control_freq;
t            = 0;

p.dt = dt; % to use in updaterobotstate
rs2DSLIP = [];  % will be initialized after first LCM message

while true
    msg = aggregator.getNextMessage(0);
    if isempty(msg)
        rate_ctrl.waitfor();
        t = t + dt;
        continue;
    end

    lc_state = eval("lcm_msgs."+lcm_state_topic+"_t(msg.data)");
    lc_cmd   = eval("lcm_msgs."+lcm_cmd_topic+"_t()");

    % update full robot state, rs = 13x1
    rs = updateRobotState(t, lc_state, rs, p);

    % Initialize SLIP leg state on first iteration
    if isempty(rs2DSLIP)
        rs2DSLIP = init2DSLIPState(rs, p);
    end

    % update individual leg SLIP dynamics (phase and forces)
    [rs2DSLIP, F_legs] = update2DSLIPState(rs2DSLIP, rs, p, dt, t);

    % stance/flight scheduling from individual leg SLIP states
    rs.ifStance = [double(rs2DSLIP.legs(1).inStance); double(rs2DSLIP.legs(2).inStance)];

    % get desired state from SLIP trajectory
    rsDes = updateRobotStateDes_longjump(rs, p, t, jump_traj);
    
    % Add SLIP leg forces to desired state for WBC to track
    rsDes.F_slip = F_legs;  % [Fx_R; Fz_R; Fx_L; Fz_L]

    % WBC: compute torques
    [tau, ddq, F] = my_WBC(controller_WBC, rs, rsDes, p); 

    % publish commands
    lc_cmd.qj_tau = tau;
    lc.publish(char(lcm_cmd_topic), lc_cmd);

    rate_ctrl.waitfor();
    t = t + dt;
end

%% functions

function jump_traj = slip_traj_opt(p, rs0)
% slip_traj_opt  direct collocation (trapezoidal) for a single long jump
%   jump_traj.t  1 x (N+1) time grid
%   jump_traj.X  4 x (N+1) states
%   jump_traj.U  2 x (N+1) controls

    import casadi.*
    N  = p.N;
    tf = p.tf;
    dt = tf / N;
    X0 = rs0.X0; % init state

    opti = casadi.Opti();

    X = opti.variable(4, N+1);   % states
    U = opti.variable(2, N+1);   % controls

    % --- initial condition ---
    opti.subject_to(X(:,1) == X0);

    obj = MX(0);

    for kk = 1:N
        Xk   = X(:,kk);
        Uk   = U(:,kk);
        Xkp1 = X(:,kk+1);
        Ukp1 = U(:,kk+1);

        % ----- dynamics (trapezoidal collocation) -----
        fk   = dynamics_slip_cart(Xk,   Uk,   p);
        fkp1 = dynamics_slip_cart(Xkp1, Ukp1, p);

        opti.subject_to(Xkp1 == Xk + dt/2 * (fk + fkp1));

        % ----- path constraints -----
        sigma = p.contact_mask(kk);   % 1=stance, 0=flight

        Fx_max = p.Umax(1);
        Fz_max = p.Umax(2);

        % vertical GRF (only in stance)
        opti.subject_to(0 <= Uk(2) <= sigma * Fz_max);

        % horizontal GRF (only in stance)
        opti.subject_to(-sigma * Fx_max <= Uk(1) <= sigma * Fx_max);

        % CoM above ground
        opti.subject_to(Xk(2) >= p.z_min);

        % ----- cost: effort + smoothness -----
        obj = obj + 1e-3 * (Uk.' * Uk);                    % effort
        obj = obj + 1e-3 * (Ukp1 - Uk).' * (Ukp1 - Uk);    % smoothness
    end

    % --- terminal constraints for one jump ---
    z0 = X0(2);

    x_final  = X(1,end);
    dx_final = X(3,end);
    dz_final = X(4,end);

    % land back at similar CoM height
    opti.subject_to(X(2,end) >= 0.9*z0);
    opti.subject_to(X(4,end) == 0);

    % --- terminal cost: maximize distance, penalize final velocity ---
    obj = obj - p.w_dist * x_final + p.w_vf * (dx_final^2 + dz_final^2);

    opti.minimize(obj);

    % --- initial guesses ---
    opti.set_initial(X, repmat(X0, 1, N+1));
    opti.set_initial(U, 0);

    opti.solver('ipopt', struct(), struct('print_level', 0));
    sol = opti.solve();

    jump_traj.t = linspace(0, tf, N+1);
    jump_traj.X = sol.value(X);
    jump_traj.U = sol.value(U);
end

function dX = dynamics_slip_cart(X, U, p)
    x   = X(1); 
    z   = X(2); 
    dx  = X(3);
    dz  = X(4);

    Fx = U(1);
    Fz = U(2);

    m = p.M_eff;
    g = p.g;

    ddx = Fx / m;
    ddz = (Fz - m*g) / m;

    dX = [dx; dz; ddx; ddz];
end

function rs0 = initRobotState(p)
    rs0.ptor = [0; p.z_min; 0];   % [x; z; pitch]
    rs0.vtor = [2.0; 2.0; 0];  % start with initial forward/upward velocity
    rs0.X0 = [0; 2.0; p.z_min; 2.0]; 
    rs0.phase = [1; 1];      % both legs start in stance
end

function rs = updateRobotState(t, state, rs, p)

    pos  = state.position;
    vel  = state.velocity;
    rpy  = state.rpy;
    omega = state.omega;
    qj   = state.qj_pos;
    dqj  = state.qj_vel;

    rs.qb  = [pos(1); pos(3); rpy(2)];        % [x; z; pitch]
    rs.dqb = [vel(1); vel(3); omega(2)];
    rs.acc = state.acceleration([1,3]);
    rs.qj  = qj;
    rs.dqj = dqj;
    rs.q   = [rs.qb;  rs.qj];
    rs.dq  = [rs.dqb; rs.dqj];
    
    % Force hip roll joints to exactly 0 to avoid NaN in dynamics
    % q = [x, z, pitch, left_hip_yaw, left_hip_roll, left_hip_pitch, left_knee, left_ankle,
    %      right_hip_yaw, right_hip_roll, right_hip_pitch, right_knee, right_ankle]
    % rs.q(5) = 0;   % left_hip_roll is q5 (4th joint after base)
    % rs.q(10) = 0;  % right_hip_roll is q10 (9th joint after base)
    % rs.dq(5) = 0;
    % rs.dq(10) = 0;

    [ptor, Jtor, dJtordq] = fcn_torso_p_J_dJdq(rs.q, rs.dq, p.l);
    [pf,   Jf,   dJfdq]   = fcn_foot_p_J_dJdq(rs.q, rs.dq, p.l);
    [H, bias]             = fcn_droid_Mass_bias(rs.q, rs.dq, p.params);

    rs.ptor   = ptor;
    rs.Jtor   = Jtor;
    rs.vtor   = Jtor * rs.dq;
    rs.dJtordq = dJtordq;

    rs.pf   = pf;
    rs.Jf   = Jf;
    rs.vf   = Jf * rs.dq;
    rs.dJfdq = dJfdq;

    rs.H    = H;
    rs.bias = bias;

    % relative torso->foot vectors (optional, obs)
    rc2f = [pf(1:2)-ptor(1:2), pf(3:4)-ptor(1:2)];
    vc2f = [rs.vf(1:2)-rs.vtor(1:2), rs.vf(3:4)-rs.vtor(1:2)];
    rs.obs = [-rc2f(2,:); -vc2f];

    % --- gait / stance logic ---

    % Initialize phases once
    if t < p.dt
        rs.phase    = [0; 0.5];                 % out-of-phase legs
        rs.pf_trans = reshape(rs.pf, [2,2]);
    end

    % If using SLIP-based long jump, we override ifStance in main loop.
    if isfield(p, 'use_slip') && p.use_slip
        % Do not update rs.ifStance here; main loop will set it
        return;
    end

    % Otherwise, original walking phase machine
    phase  = rs.phase;
    dphase = p.dt / p.T;

    for i_leg = 1:2
        idx = (i_leg-1) * 2 + (1:2);
        if phase(i_leg) < 0.5 && phase(i_leg) + dphase > 0.5
            rs.pf_trans(:,i_leg) = rs.pf(idx);
        end
    end

    rs.phase    = mod(phase + dphase, 1);
    rs.ifStance = double(rs.phase < 0.5);
end

function rsDes = updateRobotStateDes_longjump(rs, p, t, jump_traj)
% update robot torso desired pos based on offline traj
%   rs      - current robot state (struct)
%   p       - parameter structure (for future tuning, unused here)
%   t       - current time since motion start [s]
%   jump_traj.t, .X - SLIP time grid and states [x; z; dx; dz]

    rsDes = rs;  % start from current as baseline

    t_vec = jump_traj.t;
    Xref  = jump_traj.X;   % [x; z; dx; dz] 4 x (N+1)

    x_ref  = interp1(t_vec, Xref(1,:), t, 'linear', 'extrap');
    z_ref  = interp1(t_vec, Xref(2,:), t, 'linear', 'extrap');
    dx_ref = interp1(t_vec, Xref(3,:), t, 'linear', 'extrap');
    dz_ref = interp1(t_vec, Xref(4,:), t, 'linear', 'extrap');

    % --- torso position/velocity targets (x,z) ---
    rsDes.ptor(1) = x_ref;
    rsDes.ptor(2) = z_ref;
    % keep torso pitch as current or set to nominal:
    % rsDes.ptor(3) = 0;  % e.g. upright

    rsDes.vtor(1) = dx_ref;
    rsDes.vtor(2) = dz_ref;

    % --- joint posture: keep near initial pose ---
    if ~isfield(p, 'qj0')
        p.qj0  = rs.qj;
        p.dqj0 = zeros(size(rs.dqj));
    end
    rsDes.qj  = p.qj0;
    rsDes.dqj = p.dqj0;

    % --- foot velocities: zero during stance, free during flight ---
    % let feet slide slightly as torso accelerates forward
    rsDes.pf = rs.pf;  % desired position = current (no position error)
    rsDes.vf = zeros(size(rs.vf));  % desired velocity = zero (damp motion)
    
    % During flight, release foot tasks entirely (already zeroed in my_WBC)
end

%% ===WBC

function controller_WBC = ctrl_WBC(p)

    l      = p.l;
    params = p.params;
    mu     = p.mu;
    alpha  = p.alpha; 

    N    = 13;  % generalized coords (base + 10 joints)
    Ntau = 10;  % actuated joints

    import casadi.*
    opti = casadi.Opti('conic');

    % variables
    ddq = opti.variable(N,1);
    tau = opti.variable(Ntau,1);
    F   = opti.variable(4,1);

    % parameters
    q       = opti.parameter(N,1);
    dq      = opti.parameter(N,1);
    Ntask   = 3 + 4;           % 3 torso + 4 feet tasks
    Atask   = opti.parameter(Ntask,N);
    btask   = opti.parameter(Ntask,1);
    contact = opti.parameter(2,1);

    input = {q, dq, Atask, btask, contact};

    [pf, Jf, dJfdq] = fcn_foot_p_J_dJdq(q, dq, l); 
    [H, bias]       = fcn_droid_Mass_bias(q, dq, params);
    B               = [zeros(3,Ntau); eye(Ntau)];

    % objective
    obj   = MX(0);
    e_tor = Atask * ddq + btask;
    obj   = obj + e_tor.' * e_tor;
    obj   = obj + ddq.' * 1e-5 * ddq;
    obj   = obj + F.'   * 1e-4 * F;

    opti.minimize(obj);

    % dynamics constraints
    opti.subject_to(H * ddq + bias == B * tau + Jf.' * F);

    % contact-scaled friction cone and normal force bounds
    F_R = F(1:2);
    F_L = F(3:4);
    % enforce nonnegative normal forces in stance, zero in flight
    opti.subject_to(0 <= F_R(2) <= contact(1) * 500);
    opti.subject_to(0 <= F_L(2) <= contact(2) * 500);
    % friction cone ties Fx to Fz (Fx -> 0 when Fz -> 0 in flight)
    opti.subject_to(-mu * F_R(2) <= F_R(1) <= mu * F_R(2));
    opti.subject_to(-mu * F_L(2) <= F_L(1) <= mu * F_L(2));

    output = {tau, ddq, F};

    opti.solver('osqp');

    controller_WBC = opti.to_function('WBC', input, output);
end

function [tau, ddq, F] = my_WBC(ctrl_WBC, rs, rsDes, p)

    % torso task
    Kpt = diag([0 100 80]);  % increased x gain for forward tracking
    Kdt = diag([5 5 10]);

    dvTorDes = Kpt * (rsDes.ptor - rs.ptor) + ...
               Kdt * (rsDes.vtor - rs.vtor);

    Atask = rs.Jtor;
    Atask = diag(p.Q_tor) * Atask;
    btask = rs.dJtordq - dvTorDes;
    btask = diag(p.Q_tor) * btask;

    % foot
    Kpf = diag(repmat([500 800], [1,2]));
    Kdf = diag(repmat([5 5], [1,2]));

    dvFootDes = Kpf * (rsDes.pf - rs.pf) + Kdf * (rsDes.vf - rs.vf);

    Atask = [Atask; diag(p.Q_f) * rs.Jf];
    btask = [btask; diag(p.Q_f) * (rs.dJfdq - dvFootDes)];

    % joint, adjusted for 5 dof legs 
    Kp_j = diag(5 * [1 1 1 1 1 1 1 1 1 1]); % 10 dofs in total for leg
    Kd_j = diag(1 * [1 1 1 1 1 1 1 1 1 1]);

    ddqDes = Kp_j * (rsDes.qj - rs.qj) + Kd_j * (rsDes.dqj - rs.dqj);

    Q_ddq = diag(0.1 * [1 1 1 1 1 1 1 1 1 1]);
    for i_leg = 1:2
        idx = (i_leg-1) * 5 + (1:5);
        if rs.ifStance(i_leg)
            Q_ddq(idx,idx) = zeros(5);
        end
    end

    [tau, ddq, F] = ctrl_WBC(rs.q, rs.dq, Atask, btask, rs.ifStance);

    tau = full(tau);
    ddq = full(ddq);
    F   = full(F);
end

function out = hatMap(in)
    out = [0 -in(3) in(2);
           in(3) 0 -in(1);
           -in(2) in(1) 0];
end

function rs2DSLIP = init2DSLIPState(rs, p)
% init state for 2 slip legs 
%   each leg has: length, angle, velocity, stance flag, foot position

    % COM position from torso (assuming torso ≈ COM)
    rs2DSLIP.x = rs.ptor(1);
    rs2DSLIP.z = rs.ptor(2);
    rs2DSLIP.dx = rs.vtor(1);
    rs2DSLIP.dz = rs.vtor(2);
    
    % Initialize each leg
    for i = 1:2
        % Foot position: [right_x, right_z, left_x, left_z]
        foot_idx = (i-1)*2 + (1:2);
        foot_pos = rs.pf(foot_idx);
        
        % Leg vector: from foot to COM
        dx_leg = rs2DSLIP.x - foot_pos(1);
        dz_leg = rs2DSLIP.z - foot_pos(2);
        
        rs2DSLIP.legs(i).length = sqrt(dx_leg^2 + dz_leg^2);
        rs2DSLIP.legs(i).angle = atan2(dz_leg, dx_leg);  % angle from foot to COM
        rs2DSLIP.legs(i).dlength = 0;
        rs2DSLIP.legs(i).dangle = 0;
        rs2DSLIP.legs(i).footPos = foot_pos;
        
        % State transition flags
        rs2DSLIP.legs(i).inStance = true;  % start in stance
        rs2DSLIP.legs(i).prevLength = rs2DSLIP.legs(i).length;
        rs2DSLIP.legs(i).stance_start_time = 0;
        
        % Touchdown angle (for foot placement)
        rs2DSLIP.legs(i).touchdown_angle = rs2DSLIP.legs(i).angle;
    end
    
    % Global parameters
    rs2DSLIP.g = p.g;
    rs2DSLIP.M = p.M_eff;  % total mass
end

function [rs2DSLIP, F_legs] = update2DSLIPState(rs2DSLIP, rs, p, dt, t)
% Update individual leg SLIP dynamics and compute leg forces
%   Inputs:
%       rs2DSLIP - SLIP leg state
%       rs       - full robot state (for current COM/foot positions)
%       p        - parameters
%       dt       - timestep
%       t        - current time
%   Outputs:
%       rs2DSLIP - updated SLIP state
%       F_legs   - [Fx_R; Fz_R; Fx_L; Fz_L] leg forces

    % Update COM from robot state
    rs2DSLIP.x = rs.ptor(1);
    rs2DSLIP.z = rs.ptor(2);
    rs2DSLIP.dx = rs.vtor(1);
    rs2DSLIP.dz = rs.vtor(2);
    
    Fx_tot = 0;
    Fz_tot = 0;
    F_legs = zeros(4,1);
    
    for i = 1:2
        % Get current foot position from robot state
        foot_idx = (i-1)*2 + (1:2);
        current_foot_pos = rs.pf(foot_idx);
        
        % Update leg kinematics
        dx_leg = rs2DSLIP.x - current_foot_pos(1);
        dz_leg = rs2DSLIP.z - current_foot_pos(2);
        
        new_length = sqrt(dx_leg^2 + dz_leg^2);
        new_angle = atan2(dz_leg, dx_leg);
        
        % Leg velocity
        rs2DSLIP.legs(i).dlength = (new_length - rs2DSLIP.legs(i).prevLength) / dt;
        rs2DSLIP.legs(i).length = new_length;
        rs2DSLIP.legs(i).angle = new_angle;
        rs2DSLIP.legs(i).prevLength = new_length;
        
        % --- SLIP State Transitions ---
        
        % Flight -> Stance transition
        if (~rs2DSLIP.legs(i).inStance) && ...
           (rs2DSLIP.legs(i).length <= p.l0) && ...
           (rs2DSLIP.dz < 0)  % descending
            
            rs2DSLIP.legs(i).inStance = true;
            rs2DSLIP.legs(i).stance_start_time = t;
            
            % Lock foot position at touchdown
            rs2DSLIP.legs(i).footPos = current_foot_pos;
            rs2DSLIP.legs(i).touchdown_angle = new_angle;
        end
        
        % --- Stance phase: compute leg forces ---
        Fx_leg = 0;
        Fz_leg = 0;
        
        if rs2DSLIP.legs(i).inStance
            % Spring force (compression is positive force)
            F_spring = p.k * (p.l0 - rs2DSLIP.legs(i).length);
            
            % Damping force (opposes compression velocity)
            F_damp = -p.c * rs2DSLIP.legs(i).dlength;
            
            % Total leg force magnitude (along leg axis)
            F_leg_mag = F_spring + F_damp;
            
            % Only apply positive (pushing) forces
            if F_leg_mag > 0
                Fx_leg = F_leg_mag * cos(rs2DSLIP.legs(i).angle);
                Fz_leg = F_leg_mag * sin(rs2DSLIP.legs(i).angle);
            end
            
            % Stance -> Flight transition
            if (rs2DSLIP.legs(i).length >= p.l0) && ...
               (rs2DSLIP.legs(i).dlength > 0)  % extending
                
                rs2DSLIP.legs(i).inStance = false;
            end
        end
        
        % Accumulate forces
        Fx_tot = Fx_tot + Fx_leg;
        Fz_tot = Fz_tot + Fz_leg;
        
        % Store individual leg forces [Right, Left]
        if i == 1  % Right leg
            F_legs(1) = Fx_leg;
            F_legs(2) = Fz_leg;
        else  % Left leg
            F_legs(3) = Fx_leg;
            F_legs(4) = Fz_leg;
        end
    end
    
end

function FK_FD_droid()
    syms q [13,1] real
    syms dq [13,1] real
    
    syms g [1,1] real
    
    syms l_hip_roll l_hip_pitch l_thigh l_shin l_foot real
    syms M_trunk M_thigh M_shin M_foot real % Mass
    syms I_trunk I_thigh I_shin I_foot real % Moment of Inertia - resistance to rotation
    
    l = [l_hip_roll; l_hip_pitch; l_thigh; l_shin; l_foot];
    M = [M_trunk; M_thigh; M_shin; M_foot]; 
    I = [I_trunk; I_thigh; I_shin; I_foot]; 
    
    params = [g;l;M;I];
    
    N = 13;
    
    rot = @(x)[cos(x) sin(x); -sin(x) cos(x)];
    
    T01 = [rot(q3) [q1; q2]; 0 0 1];
    
    % Left Leg (hmt are always 3x3 - we are in planar motion so only x z)
    T12L = [eye(2) [0; 0]; 0 0 1];
    
    T23L = [eye(2) [0; -l_hip_roll]; 0 0 1];
    
    T34L = [rot(q6) [0; 0]; 0 0 1];
    
    T45L = [eye(2) [0; -l_thigh]; 0 0 1];
    
    T56L = [rot(q7) [0; 0]; 0 0 1];
    
    T67L = [eye(2) [0; -l_shin]; 0 0 1];
    
    T78L = [rot(q8) [0; 0]; 0 0 1];
    
    T89L = [eye(2) [0; -l_foot]; 0 0 1];
    
    T09L = T01 * T12L * T23L * T34L * T45L * T56L * T67L * T78L * T89L;
    pf_L = T09L(1:2, 3);
    
    
    % Right Leg
    T12R = [eye(2) [0; 0]; 0 0 1];
    
    T23R = [eye(2) [0; -l_hip_roll]; 0 0 1];
    
    T34R = [rot(q11) [0; 0]; 0 0 1];
    
    T45R = [eye(2) [0; -l_thigh]; 0 0 1];
    
    T56R = [rot(q12) [0; 0]; 0 0 1];
    
    T67R = [eye(2) [0; -l_shin]; 0 0 1];
    
    T78R = [rot(q13) [0; 0]; 0 0 1];
    
    T89R = [eye(2) [0; -l_foot]; 0 0 1];
    
    T09R = T01 * T12R * T23R * T34R * T45R * T56R * T67R * T78R * T89R;
    pf_R = T09R(1:2, 3);
    
    pf = [pf_R; pf_L];
    Jf = jacobian(pf, q);
    dJfdq = reshape(jacobian(Jf(:),q) * dq, size(Jf)) * dq;
    
    matlabFunction(pf, Jf, dJfdq, "File","fcn_foot_p_J_dJdq.m", "Vars",{q, dq, l})
    
    % Torso kinematics
    ptor = [T01(1:2,3); q3];
    Jtor = jacobian(ptor, q);
    dJtordq = reshape(jacobian(Jtor(:),q) * dq, size(Jtor)) * dq;
    
    matlabFunction(ptor, Jtor, dJtordq, "File","fcn_torso_p_J_dJdq.m", "Vars",{q, dq, l})
    
    p_trunk = T01(1:2,3); % floating base position is defined as trunk (entire rigid upper body) com (can add offset l_trunk/2)
    
    % Left leg - COM (Center of Mass) (uniform mass distribution)
    T0thighmidL = T01 * T12L * T23L * T34L * T45L;
    pthighL = T0thighmidL(1:2,3) - [0; l_thigh/2];
    
    T0shinmidL = T0thighmidL * T56L * T67L;
    pshinmidL = T0shinmidL(1:2,3) - [0; l_shin/2];
    
    T0footmidL = T0shinmidL * T78L * T89L;
    pfootmidL = T0footmidL(1:2,3) - [0; l_foot/2];
    
    
    % Right leg - COM (Center of Mass)so
    T0thighmidR = T01 * T12R * T23R * T34R * T45R;
    pthighR = T0thighmidR(1:2,3) - [0; l_thigh/2];
    
    T0shinmidR = T0thighmidR * T56R * T67R;
    pshinmidR = T0shinmidR(1:2,3) - [0; l_shin/2];
    
    T0footmidR = T0shinmidR * T78R * T89R;
    pfootmidR = T0footmidR(1:2,3) - [0;l_foot/2];
    
    pcom = [p_trunk, pthighL, pshinmidL, pfootmidL, pthighR, pshinmidR, pfootmidR];
    pcom = simplify(pcom);
    vcom = reshape(jacobian(pcom(:),q) * dq, size(pcom));
    
    %% Mass and Inertia matrix
    M_ = [M_trunk M_thigh M_shin M_foot M_thigh M_shin M_foot]';
    I_ = [I_trunk I_thigh I_shin I_foot I_thigh I_shin I_foot]';
    PE = pcom(2,:) * M_ * g;
    
    w_ = [dq3 dq3+dq6+dq7+dq8 dq3+dq6+dq7+dq8 dq3+dq11+dq12+dq13 dq3+dq11+dq12+dq13]';
    
    KE = sym(0);
    for ii = 1:5
        KE = KE + 0.5 * vcom(:,ii)' * M_(ii) * vcom(:,ii);
        KE = KE + 0.5 * w_(ii)' * I_(ii) * w_(ii);
    end
    
    G = jacobian(PE, q).';
    H = simplify(hessian(KE,dq));
    
    syms C [N, N] real
    for k = 1:N
        for j = 1:N
            C(k,j) = 0;
            for i = 1:N
                xkji = jacobian(H(k,j),q(i));
                xkij = jacobian(H(k,i),q(j));
                xijk = jacobian(H(i,j),q(k));
                C(k,j) = C(k,j) + 1/2 * (xkji+xkij-xijk) * dq(i);
            end
        end
    end
    
    bias = C * dq + G;
    
    matlabFunction(H, bias, 'File', 'fcn_droid_Mass_bias', 'Vars',{q,dq,params})

end