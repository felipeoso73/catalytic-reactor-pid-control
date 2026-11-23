%% =========================================================================
%  PID TEMPERATURE CONTROL IN A TWO-PHASE CATALYTIC REACTOR (CSTR)
%  Process: Catalytic valorization of agro-industrial waste at P = 25 bar
%  Nominal conversion: 73.0% at T_sp = 415.0 K
%  Primary instrumentation: RTD Pt100 Class A sensor (4-wire, 316L thermowell)
%  Compatibility: MATLAB and GNU Octave
%  =========================================================================
clear; clc; close all;

%% 1. PHYSICO-CHEMICAL PARAMETERS AND NOMINAL OPERATING CONDITIONS
P_op     = 25.0;        % Constant two-phase system pressure [bar]
V        = 2.5;         % Effective liquid-catalyst slurry volume [m^3]
F        = 0.005;       % Volumetric feed flow rate [m^3/s]
tau      = V / F;       % Nominal space time [s] (500 s)

rho      = 850.0;       % Average reactive mixture density [kg/m^3]
Cp       = 2800.0;      % Specific heat capacity [J/(kg*K)]
dHr      = -95000.0;    % Exothermic heat of reaction [J/mol]
UA       = 1800.0;      % Overall heat transfer coefficient-area product [W/K]

Ea       = 65000.0;     % Apparent activation energy (Ni-Mo/Al2O3) [J/mol]
R        = 8.314;       % Universal gas constant [J/(mol*K)]

T_sp     = 415.0;       % Reactor temperature setpoint [K]
X_nom    = 0.73;        % Target nominal conversion efficiency [-] (73.0%)

% Exact calibration of pre-exponential factor k0 for X_nom = 73.0% at 415 K
k_nom    = (X_nom / (1.0 - X_nom)) / tau;
k0       = k_nom * exp(Ea / (R * T_sp));

% Nominal feed and steady-state operating conditions
CA0_nom  = 1200.0;      % Nominal waste substrate feed concentration [mol/m^3]
T0_nom   = 385.0;       % Nominal feed temperature [K]
CA_nom   = CA0_nom * (1.0 - X_nom);

% Steady-state energy balance to compute nominal jacket temperature Tc
Q_conv   = rho * Cp * F * (T0_nom - T_sp);
Q_rxn    = (-dHr) * V * k_nom * CA_nom;
Tc_nom   = T_sp - (Q_conv + Q_rxn) / UA;

%% 2. PRIMARY INSTRUMENTATION AND PID CONTROLLER PARAMETERS
% Industrial RTD Pt100 Class A sensor inside AISI 316L tapered thermowell
tau_sens = 8.0;         % Thermal time constant of sensor assembly [s]

% PID controller tuning (reverse-acting control on jacket temperature Tc)
Kc       = 4.5;         % Proportional gain [K/K]
tau_I    = 140.0;       % Integral time constant [s]
tau_D    = 12.0;        % Derivative time constant [s]

% Physical saturation limits of the jacket heat-transfer fluid
Tc_min   = 320.0;       % Minimum cooling temperature [K]
Tc_max   = 430.0;       % Maximum heating temperature [K]

%% 3. DYNAMIC SIMULATION SETUP AND LOAD DISTURBANCES
dt       = 1.0;         % Sampling and integration time step [s]
t_end    = 3600.0;      % Simulation horizon [s] (60 minutes)
t        = 0:dt:t_end;
N        = length(t);

% Storage arrays: Open-Loop (OL) vs. Closed-Loop PID (CL)
T_ol     = zeros(1, N);
X_ol     = zeros(1, N);

T_cl     = zeros(1, N);
Tm_cl    = zeros(1, N);
X_cl     = zeros(1, N);
Tc_cl    = zeros(1, N);

% Initial steady-state vector: [CA; T; Tm]
state_ol = [CA_nom; T_sp; T_sp];
state_cl = [CA_nom; T_sp; T_sp];

integral_err = 0.0;
prev_Tm      = T_sp;

% Coupled differential equations (mass balance, energy balance, sensor lag)
ode_rhs = @(s, CA0, T0, Tc) [ ...
    (F / V) * (CA0 - s(1)) - (k0 * exp(-Ea / (R * s(2)))) * s(1); ...
    (F / V) * (T0 - s(2)) + ((-dHr) / (rho * Cp)) * (k0 * exp(-Ea / (R * s(2)))) * s(1) ...
        - (UA / (rho * Cp * V)) * (s(2) - Tc); ...
    (s(2) - s(3)) / tau_sens ...
];

%% 4. NUMERICAL INTEGRATION LOOP (RK4) AND DISCRETE PID CONTROL
for i = 1:N
    ti = t(i);
    
    % Disturbance 1 (t = 600 s): +20% step increase in waste feed concentration
    if ti >= 600
        CA0 = 1.20 * CA0_nom;
    else
        CA0 = CA0_nom;
    end
    
    % Disturbance 2 (t = 1800 s): -10 K step drop in feed temperature
    if ti >= 1800
        T0 = T0_nom - 10.0;
    else
        T0 = T0_nom;
    end
    
    % --- Open-Loop Simulation (Tc fixed at Tc_nom) ---
    T_ol(i) = state_ol(2);
    X_ol(i) = ((CA0 - state_ol(1)) / CA0) * 100.0;
    
    k1 = ode_rhs(state_ol, CA0, T0, Tc_nom);
    k2 = ode_rhs(state_ol + 0.5*dt*k1, CA0, T0, Tc_nom);
    k3 = ode_rhs(state_ol + 0.5*dt*k2, CA0, T0, Tc_nom);
    k4 = ode_rhs(state_ol + dt*k3, CA0, T0, Tc_nom);
    state_ol = state_ol + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
    
    % --- Closed-Loop Simulation (PID Control with RTD Pt100) ---
    T_cl(i)  = state_cl(2);
    Tm_cl(i) = state_cl(3);
    X_cl(i)  = ((CA0 - state_cl(1)) / CA0) * 100.0;
    
    err      = T_sp - Tm_cl(i);
    dTm_dt   = (Tm_cl(i) - prev_Tm) / dt;
    prev_Tm  = Tm_cl(i);
    
    % PID control law (derivative-on-measurement to prevent derivative kick)
    u_unsat  = Tc_nom + Kc * (err + (1.0 / tau_I) * integral_err - tau_D * dTm_dt);
    Tc_cmd   = min(max(u_unsat, Tc_min), Tc_max);
    
    % Conditional integration (Anti-Windup)
    if (u_unsat == Tc_cmd) || ((u_unsat > Tc_max) && (err < 0)) || ((u_unsat < Tc_min) && (err > 0))
        integral_err = integral_err + err * dt;
    end
    
    Tc_cl(i) = Tc_cmd;
    
    k1 = ode_rhs(state_cl, CA0, T0, Tc_cmd);
    k2 = ode_rhs(state_cl + 0.5*dt*k1, CA0, T0, Tc_cmd);
    k3 = ode_rhs(state_cl + 0.5*dt*k2, CA0, T0, Tc_cmd);
    k4 = ode_rhs(state_cl + dt*k3, CA0, T0, Tc_cmd);
    state_cl = state_cl + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
end

%% 5. PERFORMANCE METRICS AND CONSOLE SUMMARY
IAE_ol = trapz(t, abs(T_sp - T_ol));
IAE_cl = trapz(t, abs(T_sp - T_cl));

fprintf('=================================================================\n');
fprintf(' PERFORMANCE SUMMARY - TWO-PHASE CATALYTIC REACTOR (P = %.1f bar)\n', P_op);
fprintf('=================================================================\n');
fprintf(' Initial nominal conversion        : %.2f %%\n', X_cl(1));
fprintf(' Peak temperature (Open-Loop)      : %.2f K\n', max(T_ol));
fprintf(' Peak temperature (PID Control)    : %.2f K\n', max(T_cl));
fprintf(' Final conversion (Open-Loop)      : %.2f %%\n', X_ol(end));
fprintf(' Final conversion (PID Control)    : %.2f %%\n', X_cl(end));
fprintf(' Integral Absolute Error IAE (OL)  : %.1f K*s\n', IAE_ol);
fprintf(' Integral Absolute Error IAE (PID) : %.1f K*s (Reduction: %.1f%%)\n', ...
        IAE_cl, (1 - IAE_cl/IAE_ol)*100);
fprintf('=================================================================\n');

%% 6. COMPARATIVE PLOTS AND HIGH-RESOLUTION EXPORT FOR PORTFOLIO
t_min = t / 60.0;
fig = figure('Color', 'w', 'Position', [100, 100, 950, 720]);

% Subplot 1: Reactor temperature
subplot(3, 1, 1);
plot(t_min, T_ol, '--r', 'LineWidth', 1.6); hold on;
plot(t_min, T_cl, '-b',  'LineWidth', 1.8);
plot(t_min, ones(size(t_min))*T_sp, ':k', 'LineWidth', 1.2);
ylabel('Temperature [K]', 'FontWeight', 'bold');
title('PID Temperature Control in a Two-Phase Catalytic Reactor (P = 25 bar)', 'FontSize', 11);
legend('Open-Loop', 'Closed-Loop (PID + Pt100)', 'Setpoint (415 K)', 'Location', 'best');
grid on; xlim([0 60]);

% Subplot 2: Conversion efficiency
subplot(3, 1, 2);
plot(t_min, X_ol, '--r', 'LineWidth', 1.6); hold on;
plot(t_min, X_cl, '-g',  'LineWidth', 1.8);
plot(t_min, ones(size(t_min))*73.0, ':k', 'LineWidth', 1.2);
ylabel('Conversion X_A [%]', 'FontWeight', 'bold');
legend('Open-Loop', 'Closed-Loop (PID)', 'Nominal Target (73.0%)', 'Location', 'best');
grid on; xlim([0 60]);

% Subplot 3: Manipulated variable (Jacket temperature)
subplot(3, 1, 3);
plot(t_min, Tc_cl, '-m', 'LineWidth', 1.6); hold on;
plot(t_min, ones(size(t_min))*Tc_nom, ':k', 'LineWidth', 1.2);
xlabel('Operating time [min]', 'FontWeight', 'bold');
ylabel('Jacket Temp. T_c [K]', 'FontWeight', 'bold');
legend('Manipulated Signal T_c (PID)', 'Nominal T_c', 'Location', 'best');
grid on; xlim([0 60]);

% Export high-resolution PNG figure for README.md and Notion
print(fig, 'reactor_pid_results_en.png', '-dpng', '-r300');
