%% v3
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

% parameters:
params.M = 80; % mass (kg)
params.g = [0; 0; -9.81];
params.l0 = 1.0;            % rest spring leg length, at TD l0 = lh
params.lh = 1.0;            % get from xml - humanoid virtual leg length used to map to SLIP leg
params.yhip = 0.0;          % 0 for now for testing - hip offset in y-dir. nominal val is torso width/2
params.th0 = deg2rad(24);   % init TD angle guess
params.ks0 = 6000;
params.tf = 5;              % single step time interval
% potential param to add: scaling param for phi (affects sagittal dir)

% %% Main simulation loop
% initial apex state [h, vx, vy]
X0 = [2.0; 3.5; 0];

% optimize for periodic gait timings
fprintf('Finding periodic gait...\n');
[X0_star, u0_star] = find_periodic_gait(X0, params);
fprintf('Periodic gait found:\n');
fprintf('  Apex height h0     = %.4f m\n', X0_star(1));
fprintf('  Forward velocity vx = %.4f m/s\n', X0_star(2));
fprintf('  Lateral velocity vy = %.4f m/s\n', X0_star(3));
fprintf('  θ   = %.3f deg\n', rad2deg(u0_star(1)));
fprintf('  φ     = %.3f deg\n', rad2deg(u0_star(2)));
fprintf('  ks1 = %.2f N/m\n', u0_star(3));
fprintf('  ks2 = %.2f N/m\n', u0_star(4));
fprintf('--------------------------------------------------\n');

% simulate
N = 5; % number of steps
X0 = X0_star;
traj = zeros(3, N);
fprintf('Starting simulation...\n');
for n = 1:N
    K = compute_deadbeat(X0_star, u0_star, params);
    u = u0_star + K * (X0 - X0_star);
    [X1, t_TD, t_LO] = slip_return_map(X0, u, params);
    traj(:, n) = X0;
    X0 = X1;
end
save('SLIP3D_data.mat', 'traj', 'K', 'X0_star', 'u0_star', 'params');
fprintf('Simulation complete.\n');

% %% Functions
function pf = get_TD_pos(X, u, params)
% eq. 3
% return: foot pos as TD happens
% input: X = [h; vx; vy]
    h = X(1);
    theta = u(1);
    phi = u(2);
    ps = [0; 0; h];                 % position of mass 3x1
    phip = [0; params.yhip; 0];     % position of hip wrt CoM = offset in y-dir 3x1
    th = u(1);
    phi = u(2); 
    lh = params.lh;

    pf = ps + phip + lh * [sin(th)*cos(phi); -sin(th)*sin(phi); -cos(th)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = dynamics_SLIP(~, X, phase, params, pf, ks) % time not relevant
% pf = foot pos at most recent TD (to allow for both flight/stance scenarios)
% ks = will be ks1 before max compression, ks2 after max compression
    M = params.M;
    g = params.g; 
    l0 = params.l0;
    p = X(1:3);
    pd = X(4:6);
    if strcmp(phase,'flight')   % ballistic dynamics
        dX = [pd; g];
    else                        % stance dynamics: eq. 2
        l = p - pf;
        lhat = l / norm(l); % use hat map or l / norm(l)?
        F = ks * (l0 - norm(l)) * lhat + M * g;
        pdd = F / M;
        dX = [pd; pdd];
    end
end

function out = hatMap(in) % not used
    out = [0 -in(3) in(2);
           in(3) 0 -in(1);
           -in(2) in(1) 0];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value, isterminal, direction] = flight_event(~, X, u, params) 
% TD event: z pos = l_h * cos(th) (eq. 3)
    ps = X(1:3);
    l0 = params.l0;
    th = u(1);
    value = ps(3) - l0 * cos(th); 
    isterminal = 1;
    direction = -1;
end

function [value, isterminal, direction] = stance_events(~, X, params, pf)
% two events to check for during stance:
% 1. max compression: l' * v = 0 and ||l|| < l0
% 2. LO event: when leg returns to rest length (eq. 4)
    ps = X(1:3);
    vs = X(4:6);
    l = ps - pf;

    value = [l.' * vs; norm(l) - params.l0];
    isterminal = [1; 1]; % currently terminating after max compression (could change)
    direction = [0; 1];
end

% NOTE: could add apex event where vertical velocity X(6) = 0 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X1, t_TD, t_LO] = slip_return_map(X, u, params)
% 4 phase integration of 1 step from X0 to X1 (apex to apex)
% returns next apex state X1, and t_TD/t_LO for cost function
% u = [th; phi; ks1; ks2]
% ks1 = during compression
% ks2 = during extension
    h = X(1); vx = X(2); vy = X(3);
    X0 = [0; 0; h; vx; vy; 0]; % expand into full SLIP state to pass into dynamics
    ks1 = u(3);
    ks2 = u(4);
    pf = get_TD_pos(X, u, params);
    tf = params.tf;

    opts_flight = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,X) flight_event(t,X,u,params)); % while in flight, detect for TD
    opts_stance = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,X) stance_events(t,X,params,pf)); % while in stance, detect for LO
    
    % phase 1 - flight: apex to TD
    fprintf('θ = %.2f deg | cos(θ) = %.3f | l0*cos(θ) = %.3f | initial height h0 = %.3f\n', ...
        rad2deg(u(1)), cos(u(1)), params.l0*cos(u(1)), X(1)); % debug

    [t1, X1, te1, Xe1] = ode45(@(t,X) dynamics_SLIP(t,X,'flight',params,pf,ks1), [0 tf], X0, opts_flight); % ks input doenst matter
    if isempty(te1)               
        fprintf('No TD detected during first flight phase'); % debug
        Xe1 = X1(end,:)';   % safeguard: just take last value to prevent crash
        t_TD = t1(end);
    else
        t_TD = te1;
    end
    fprintf('Phase 1 integration time span: [%.6f, %.6f]\n', 0, t1(end));
    fprintf('Final leg length: %.6f m\n', norm(X1(end,1:3) - pf'));
    fprintf('Final leg length rate (l''·v): %.6f\n', dot(X1(end,1:3) - pf', X1(end,4:6)));
    fprintf('Spring stiffness ks1 = %.2f N/m\n', ks1);
    fprintf('State at end of stance:\n');
    disp(X1(end,:));

    % phase 2 - stance: TD to max compression (ks1)
    [t2, X2, te2, Xe2] = ode45(@(t,X) dynamics_SLIP(t,X,'stance',params,pf,ks1), [t1(end) tf], Xe1, opts_stance); % use t1(end) to ensure continuity
    if isempty(te2)
        fprintf('No max compression detected during stance phase\n'); % debug
        fprintf('Phase 2 integration time span: [%.6f, %.6f]\n', t1(end), tf);
        fprintf('Final leg length: %.6f m\n', norm(X2(end,1:3) - pf'));
        fprintf('Final leg length rate (l''·v): %.6f\n', dot(X2(end,1:3) - pf', X2(end,4:6)));
        fprintf('Spring stiffness ks1 = %.2f N/m\n', ks1);
        fprintf('State at end of stance:\n');
        disp(X2(end,:));
        Xe2 = X2(end,:)';
    end

    % NOTE: could consider not terminating after max compression and having this 
    % under an if conditional while tuning ks values
    % phase 3 - stance: max compression to LO (ks2)
    [t3, X3, te3, Xe3] = ode45(@(t,X) dynamics_SLIP(t,X,'stance',params,pf,ks2), [t2(end) tf], Xe2, opts_stance);
    if isempty(te3)
        fprintf('No LO detected during stance phase'); % debug
        Xe3 = X3(end,:)';
        t_LO = t3(end);
    else
        t_LO = te3;
    end
    
    % phase 4 - flight: LO to apex
    [t4, X4] = ode45(@(t,X) dynamics_SLIP(t,X,'flight',params,pf,ks2), [t3(end) tf], Xe3, opts_flight); % ks input doesnt matter
    
    % could replace this if add apex event detection function
    [~, idx] = max(X4(:,3)); % apex idx = max z pos of second flight phase
    apex = X4(idx,:);        % find xyz pos at apex idx 
    
    X1 = [apex(3); apex(4); apex(5)]; % convert back into simplified apex state
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X0_star, u0_star] = find_periodic_gait(X0, params)
% solve LS problem to get (X0*, u0*) optimal state-control pair that
% achieves desired gait timings given a desired forward vel vx
    % initial decision var guess
    z0 = [2.0; 0.0; params.ks0; params.th0];
    vx = X0(2); % this is kept constant
    fun = @(z) periodic_cost(z, vx, params); % decision variables (z): [h0, vy0, ks, th]
    % need max iterations?
    options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',2000,'TolX',1e-8);
    % consider upper and lower bounds
    sol = lsqnonlin(fun,z0,[],[],options); % sol = [h0*, vy0*, ks*, th*]
    
    X0_star = [sol(1); vx; sol(2)];
    u0_star = [sol(4); 0; sol(3); sol(3)];
end

function err = periodic_cost(z, vx, params)
% form symbolic cost function (eq. 13)
    % z = [h0, vy0, ks, th] 4x1 vector of decision vars
    h0  = z(1);
    vy0 = z(2);
    ks  = z(3);
    th  = z(4);
    X0 = [h0; vx; vy0];     % eq. 14
    u0 = [th; 0; ks; ks];   % eq. 15

    A = diag([1 1 -1]);     % eq. 8
    [X1, t_TD, t_LO] = slip_return_map(X0, u0, params); % symbolically integrate one step forward
    Tdes = get_des_gait_timings(X0); 
    Tcurr = [t_TD; t_LO];   % eq. 12
    
    err = [A*X0 - X1; Tdes - Tcurr]; % eq. 13
    err = err(:); % make sure col vec (can debug w this later)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        Jx(:,i) = (slip_return_map(Xp,u0_star,params) - slip_return_map(Xm,u0_star,params))/(2*dX);
    end
    for j = 1:4
        up = u0_star; up(j)=up(j)+du;
        um = u0_star; um(j)=um(j)-du;
        Ju(:,j) = (slip_return_map(X0_star,up,params) - slip_return_map(X0_star,um,params))/(2*du);
    end
    K = -pinv(Ju)*Jx;
end
