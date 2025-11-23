%% 11/21 version working with gait library and unconstrained optimization
% but getting stuck on N=6 steps
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
params.lh = 1.1;            % get from xml - humanoid virtual leg length used to map to SLIP leg
params.yhip = 0.1;          % 0 for now for testing - hip offset in y-dir. nominal val is torso width/2
params.th0 = deg2rad(24);   % init TD angle guess
params.ks0 = 6000;          % init stiffness guess
params.tf = 5.0;            % single step time interval
% potential param to add: scaling param for phi (affects sagittal dir)

% range of desired forward velocities
vx_range = linspace(3.5, 6.5, 31); % from paper

% gait library
X0_stars = zeros(3, length(vx_range));
u0_stars = zeros(4, length(vx_range));
K_all = cell(1,length(vx_range));

% optimize periodic gait timings across vx range
fprintf('Generating gait library...\n');
for i = 1:length(vx_range)
    vx_des = vx_range(i);
    X0 = [2.0; vx_des; 0]; % initial guess apex state [h, vx, vy], h and vy are guesses to be adjusted, vx is desired vel for entire traj
    
    fprintf('\n--- Speed = %.2f m/s ---\n', vx_des);
    [X0_star, u0_star] = find_periodic_gait(X0, params);
    fprintf('Periodic gait found for %.2f m/s:\n', vx_des);
    fprintf('  Apex height h0     = %.4f m\n', X0_star(1));
    fprintf('  θ = %.2f°, ks = %.0f N/m\n', rad2deg(u0_star(1)), u0_star(3));
    fprintf('  Lateral velocity vy = %.4f m/s\n', X0_star(3));
    fprintf('  φ     = %.3f deg\n', rad2deg(u0_star(2)));
    fprintf('--------------------------------------------------\n');

    K = compute_deadbeat(X0_star, u0_star, params)
    
    X0_stars(:,i) = X0_star;
    u0_stars(:,i) = u0_star;
    K_all{i} = K;
end
save('SLIP3D_gait_library.mat','vx_range','X0_stars','u0_stars','K_all','params');
fprintf('\nGait library generation complete.\nSaved to SLIP3D_gait_library.mat\n');

% %% Main simulation loop
% simulate
N = 5; % number of steps
X0 = X0_star;
traj = zeros(3, N);
fprintf('Starting simulation...\n');
for n = 1:N
    % retrieve closest K corresponding to vx from "library":
    fprintf('Retrieving from library...\n');
    vx_curr = X0(2);
    [~, idx] = min(abs(vx_range - vx_curr));
    X0_star = X0_stars(:,idx);
    u0_star = u0_stars(:,idx);
    K = K_all{idx};

    u = u0_star + K * (X0 - X0_star);                  % eq 19
    [X1, t_TD, t_LO] = slip_return_map(X0, u, params);  % integrate one step forward with adjusted control
    traj(:, n) = X0;
    X0 = X1;
end
save('SLIP3D_data.mat', 'traj', 'K', 'X0_star', 'u0_star', 'params');
fprintf('Simulation complete.\n');

% %% Functions
function pf = get_TD_pos(X, u, params)
% eq. 3
% return: foot pos as TD happens
% input: X = [h; vx; vy]. NOT full slip state
    h = X(1);
    theta = u(1);
    phi = u(2);
    ps = [0; 0; h];                 % position of mass 3x1
    phip = [0; -params.yhip; 0];     % position of hip wrt CoM = offset in y-dir 3x1
    th = u(1);
    phi = u(2); 
    lh = params.lh;

    pf = ps + phip + lh * [sin(th)*cos(phi); -sin(th)*sin(phi); -cos(th)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dX = dynamics_SLIP(~, Xs, phase, params, pf, ks) % time not relevant
% Xs = full slip state
% pf = foot pos at most recent TD (to allow for both flight/stance scenarios)
% ks = will be ks1 before max compression, ks2 after max compression
    M = params.M;
    g = params.g; 
    l0 = params.l0;
    p = Xs(1:3);
    pd = Xs(4:6);
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
% event during 1st stance phase
% max compression event: l' * v = 0
    ps = Xs(1:3);
    vs = Xs(4:6);
    l = ps - pf;
    value = l.' * vs;
    isterminal = 1; % currently terminating after max compression (could change)
    direction = 0;
end

function [value, isterminal, direction] = LOevent(~, Xs, params, pf)
% event during 2nd stance phase
% LO event: ||l|| - l0 = 0, when leg returns to rest length (eq. 4)
    ps = Xs(1:3);
    l = ps - pf;
    value = norm(l) - params.l0;
    isterminal = 1; 
    direction = 1;
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
    Xs0 = [0; 0; h; vx; vy; 0]; % expand into full SLIP state to pass into event functions / dynamics
    ks1 = u(3);
    ks2 = u(4);
    pf = get_TD_pos(X, u, params);
    tf = params.tf;

    opts_flight = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) TDevent(t,Xs,u,params)); % while in flight, detect for TD
    opts_compression = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) MCevent(t,Xs,pf)); % while in stance, detect for MC
    opts_extension = odeset('RelTol',1e-6,'AbsTol',1e-8,'Events',@(t,Xs) LOevent(t,Xs,params,pf)); % while in stance, detect for LO
    
    % fprintf('θ = %.2f deg | cos(θ) = %.3f | lh*cos(θ) = %.3f | initial height h0 = %.3f\n', ...
    %     rad2deg(u(1)), cos(u(1)), params.lh*cos(u(1)), X(1)); % debug to check that apex is higher than TD expression
    % fprintf('Initial COM height h0 = %.2f\n', Xs0(3));
    % fprintf('Rest leg length lh = %.2f\n', params.lh);
    % fprintf('Expected TD height (lh*cos(th)) = %.2f\n', params.lh*cos(u(1)));
    % fprintf('Gravity vector = [%.2f %.2f %.2f]\n', params.g);

    % phase 1 - flight: apex to TD
    [t1, Xs1, te1, Xe1] = ode45(@(t,Xs) dynamics_SLIP(t,Xs,'flight',params,pf,ks1), [0 tf], Xs0, opts_flight); % ks input doenst matter
    if isempty(te1)               
        fprintf('No TD detected during first flight phase'); % debug
        Xe1 = Xs1(end,:)';   % safeguard: just take last value to prevent crash
        t_TD = t1(end);
    else
        t_TD = te1;
    end
    % fprintf('Phase 1 integration time span: [%.6f, %.6f]\n', 0, t1(end));
    % fprintf('Final leg length: %.6f m\n', norm(X1(end,1:3) - pf'));
    % fprintf('Final leg length rate (l''·v): %.6f\n', dot(X1(end,1:3) - pf', X1(end,4:6)));
    % fprintf('Spring stiffness ks1 = %.2f N/m\n', ks1);
    % fprintf('State at end of stance:\n');
    % disp(X1(end,:));

    % phase 2 - stance: TD to max compression (ks1)
    [t2, Xs2, te2, Xe2] = ode45(@(t,Xs) dynamics_SLIP(t,Xs,'stance',params,pf,ks1), [t1(end) tf], Xe1, opts_compression); % use t1(end) to ensure continuity
    if isempty(te2)
        fprintf('No max compression detected during stance phase\n'); % debug
        fprintf('Phase 2 integration time span: [%.6f, %.6f]\n', t1(end), tf);
        fprintf('Final leg length: %.6f m\n', norm(Xs2(end,1:3) - pf'));
        fprintf('Final leg length rate (l''·v): %.6f\n', dot(Xs2(end,1:3) - pf', Xs2(end,4:6)));
        fprintf('Spring stiffness ks1 = %.2f N/m\n', ks1);
        fprintf('State at end of stance:\n');
        disp(Xs2(end,:));
        Xe2 = Xs2(end,:)';
    end

    % NOTE: could consider not terminating after max compression and having this 
    % under an if conditional while tuning ks values
    % phase 3 - stance: max compression to LO (ks2)
    [t3, Xs3, te3, Xe3] = ode45(@(t,X) dynamics_SLIP(t,X,'stance',params,pf,ks2), [t2(end) tf], Xe2, opts_extension);
    if isempty(te3)
        fprintf('No LO detected during stance phase'); % debug
        leg_length_final = norm(Xs3(end,1:3)' - pf);
        fprintf('  Final leg length: %.4f m (l0=%.4f m), diff=%.4f m\n', leg_length_final, params.l0, leg_length_final - params.l0);
        Xe3 = Xs3(end,:)'
        t_LO = t3(end);
    else
        t_LO = te3;
    end
    
    % phase 4 - flight: LO to apex
    [t4, Xs4] = ode45(@(t,X) dynamics_SLIP(t,X,'flight',params,pf,ks2), [t3(end) tf], Xe3, opts_flight); % ks input doesnt matter
    
    % could replace this if add apex event detection function
    [~, idx] = max(Xs4(:,3)); % apex idx = max z pos of second flight phase
    apex = Xs4(idx,:);        % find xyz pos at apex idx 
    
    X1 = [apex(3); apex(4); apex(5)]; % convert back into simplified apex state
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X0_star, u0_star] = find_periodic_gait(X0, params)
% solve LS problem to get (X0*, u0*) optimal state-control pair that
% achieves desired gait timings given a desired forward vel vx
% X0 = [h0; vx; vy0]
% decision variables (z): [h0, vy0, ks, th]

    z0 = [X0(1); X0(3); params.ks0; params.th0]; % initial guess (apex)
    vx = X0(2); % this is kept constant
    fun = @(z) periodic_cost(z, vx, params); 
    % need max iterations?
    options = optimoptions('lsqnonlin','Display','iter','MaxFunEvals',2000,'TolX',1e-8);
    
    (* % z = [h0, vy0, ks, th]
    lb = [1.5; -1.0; 4000; deg2rad(8)];             % lower bound (adjusted for l0=1.0, lh=1.1)
    ub = [2.5; 1.0; 50000; deg2rad(35)];            % upper bound *)
    sol = lsqnonlin(fun,z0,lb,ub,options); % sol = [h0*, vy0*, ks*, th*]
    
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
    
    % Add large penalty if LO was not detected (invalid gait)
    if ~LO_detected
        err = err + 100 * ones(size(err)); % Large penalty to discourage invalid gaits
    end
    
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
