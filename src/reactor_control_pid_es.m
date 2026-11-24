%% =========================================================================
%  CONTROL PID DE TEMPERATURA EN REACTOR CATALÍTICO BIFÁSICO (CSTR)
%  Proceso: Valorización catalítica de residuos agroindustriales a P = 25 bar
%  Conversión nominal: 73.0% a T_sp = 415.0 K
%  Instrumentación primaria: Sensor RTD Pt100 Clase A (4 hilos, termopozo 316L)
%  Compatibilidad: MATLAB y GNU Octave
%  =========================================================================
clear; clc; close all;

%% 1. PARÁMETROS FÍSICO-QUÍMICOS Y CONDICIONES NOMINALES DE OPERACIÓN
P_op     = 25.0;        % Presión constante del sistema bifásico [bar]
V        = 2.5;         % Volumen efectivo de la fase líquida-catalítica [m^3]
F        = 0.005;       % Caudal volumétrico de alimentación [m^3/s]
tau      = V / F;       % Tiempo espacial nominal [s] (500 s)

rho      = 850.0;       % Densidad media de la mezcla reactiva [kg/m^3]
Cp       = 2800.0;      % Capacidad calorífica específica [J/(kg*K)]
dHr      = -95000.0;    % Entalpía de reacción exotérmica [J/mol]
UA       = 1800.0;      % Coeficiente global de transferencia de calor [W/K]

Ea       = 65000.0;     % Energía de activación aparente (Ni-Mo/Al2O3) [J/mol]
R        = 8.314;       % Constante universal de los gases [J/(mol*K)]

T_sp     = 415.0;       % Temperatura de consigna (Setpoint) [K]
X_nom    = 0.73;        % Eficiencia de conversión nominal objetivo [-] (73.0%)

% Calibración exacta del factor pre-exponencial k0 para X_nom = 73.0% a 415 K
k_nom    = (X_nom / (1.0 - X_nom)) / tau;
k0       = k_nom * exp(Ea / (R * T_sp));

% Condiciones nominales de alimentación y estado estacionario
CA0_nom  = 1200.0;      % Concentración nominal del sustrato residual [mol/m^3]
T0_nom   = 385.0;       % Temperatura nominal de alimentación [K]
CA_nom   = CA0_nom * (1.0 - X_nom);

% Balance energético estacionario para determinar Tc nominal de la chaqueta
Q_conv   = rho * Cp * F * (T0_nom - T_sp);
Q_rxn    = (-dHr) * V * k_nom * CA_nom;
Tc_nom   = T_sp - (Q_conv + Q_rxn) / UA;

%% 2. INSTRUMENTACIÓN PRIMARIA Y PARÁMETROS DEL CONTROLADOR PID
% Sensor RTD Pt100 Clase A en termopozo AISI 316L (punta reducida)
tau_sens = 8.0;         % Constante de tiempo térmica del conjunto sensor [s]

% Sintonización del controlador PID (acción inversa sobre Tc)
Kc       = 4.5;         % Ganancia proporcional [K/K]
tau_I    = 140.0;       % Tiempo integral [s]
tau_D    = 12.0;        % Tiempo derivativo [s]

% Límites físicos del fluido térmico de la chaqueta (saturación)
Tc_min   = 320.0;       % Temperatura mínima de enfriamiento [K]
Tc_max   = 430.0;       % Temperatura máxima de calentamiento [K]

%% 3. CONFIGURACIÓN DE LA SIMULACIÓN DINÁMICA Y PERTURBACIONES DE CARGA
dt       = 1.0;         % Paso de muestreo e integración [s]
t_end    = 3600.0;      % Horizonte de simulación [s] (60 minutos)
t        = 0:dt:t_end;
N        = length(t);

% Vectores de almacenamiento: Lazo Abierto (OL) vs. Lazo Cerrado PID (CL)
T_ol     = zeros(1, N);
X_ol     = zeros(1, N);

T_cl     = zeros(1, N);
Tm_cl    = zeros(1, N);
X_cl     = zeros(1, N);
Tc_cl    = zeros(1, N);

% Estados iniciales en régimen estacionario: [CA; T; Tm]
state_ol = [CA_nom; T_sp; T_sp];
state_cl = [CA_nom; T_sp; T_sp];

integral_err = 0.0;
prev_Tm      = T_sp;

% Definición de ecuaciones diferenciales acopladas (balances + sensor)
ode_rhs = @(s, CA0, T0, Tc) [ ...
    (F / V) * (CA0 - s(1)) - (k0 * exp(-Ea / (R * s(2)))) * s(1); ...
    (F / V) * (T0 - s(2)) + ((-dHr) / (rho * Cp)) * (k0 * exp(-Ea / (R * s(2)))) * s(1) ...
        - (UA / (rho * Cp * V)) * (s(2) - Tc); ...
    (s(2) - s(3)) / tau_sens ...
];

%% 4. BUCLE DE INTEGRACIÓN NUMÉRICA (RK4) Y CONTROL PID DISCRETO
for i = 1:N
    ti = t(i);
    
    % Perturbación 1 (t = 600 s): +20% en concentración de residuo alimentado
    if ti >= 600
        CA0 = 1.20 * CA0_nom;
    else
        CA0 = CA0_nom;
    end
    
    % Perturbación 2 (t = 1800 s): caída de -10 K en temperatura de entrada
    if ti >= 1800
        T0 = T0_nom - 10.0;
    else
        T0 = T0_nom;
    end
    
    % --- Simulación en Lazo Abierto (Tc fija en Tc_nom) ---
    T_ol(i) = state_ol(2);
    X_ol(i) = ((CA0 - state_ol(1)) / CA0) * 100.0;
    
    k1 = ode_rhs(state_ol, CA0, T0, Tc_nom);
    k2 = ode_rhs(state_ol + 0.5*dt*k1, CA0, T0, Tc_nom);
    k3 = ode_rhs(state_ol + 0.5*dt*k2, CA0, T0, Tc_nom);
    k4 = ode_rhs(state_ol + dt*k3, CA0, T0, Tc_nom);
    state_ol = state_ol + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
    
    % --- Simulación en Lazo Cerrado (Control PID con RTD Pt100) ---
    T_cl(i)  = state_cl(2);
    Tm_cl(i) = state_cl(3);
    X_cl(i)  = ((CA0 - state_cl(1)) / CA0) * 100.0;
    
    err      = T_sp - Tm_cl(i);
    dTm_dt   = (Tm_cl(i) - prev_Tm) / dt;
    prev_Tm  = Tm_cl(i);
    
    % Ley de control PID (derivada sobre la variable medida para evitar patada derivativa)
    u_unsat  = Tc_nom + Kc * (err + (1.0 / tau_I) * integral_err - tau_D * dTm_dt);
    Tc_cmd   = min(max(u_unsat, Tc_min), Tc_max);
    
    % Integración condicional (Anti-Windup)
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

%% 5. CÁLCULO DE INDICADORES DE DESEMPEÑO Y RESUMEN EN CONSOLA
IAE_ol = trapz(t, abs(T_sp - T_ol));
IAE_cl = trapz(t, abs(T_sp - T_cl));

fprintf('=================================================================\n');
fprintf(' RESUMEN DE DESEMPEÑO - REACTOR CATALÍTICO BIFÁSICO (P = %.1f bar)\n', P_op);
fprintf('=================================================================\n');
fprintf(' Conversión nominal inicial        : %.2f %%\n', X_cl(1));
fprintf(' Temperatura máxima (Lazo Abierto) : %.2f K\n', max(T_ol));
fprintf(' Temperatura máxima (Control PID)  : %.2f K\n', max(T_cl));
fprintf(' Conversión final (Lazo Abierto)   : %.2f %%\n', X_ol(end));
fprintf(' Conversión final (Control PID)    : %.2f %%\n', X_cl(end));
fprintf(' Error Integral Absoluto IAE (OL)  : %.1f K*s\n', IAE_ol);
fprintf(' Error Integral Absoluto IAE (PID) : %.1f K*s (Reducción: %.1f%%)\n', ...
        IAE_cl, (1 - IAE_cl/IAE_ol)*100);
fprintf('=================================================================\n');

%% 6. GENERACIÓN DE GRÁFICAS COMPARATIVAS Y EXPORTACIÓN PARA PORTAFOLIO
t_min = t / 60.0;
fig = figure('Color', 'w', 'Position', [100, 100, 950, 720]);

% Subgráfica 1: Temperatura del reactor
subplot(3, 1, 1);
plot(t_min, T_ol, '--r', 'LineWidth', 1.6); hold on;
plot(t_min, T_cl, '-b',  'LineWidth', 1.8);
plot(t_min, ones(size(t_min))*T_sp, ':k', 'LineWidth', 1.2);
ylabel('Temperatura [K]', 'FontWeight', 'bold');
title('Control PID de Temperatura en Reactor Catalítico Bifásico (P = 25 bar)', 'FontSize', 11);
legend('Lazo Abierto', 'Lazo Cerrado (PID + Pt100)', 'Setpoint (415 K)', 'Location', 'best');
grid on; xlim([0 60]);

% Subgráfica 2: Eficiencia de conversión
subplot(3, 1, 2);
plot(t_min, X_ol, '--r', 'LineWidth', 1.6); hold on;
plot(t_min, X_cl, '-g',  'LineWidth', 1.8);
plot(t_min, ones(size(t_min))*73.0, ':k', 'LineWidth', 1.2);
ylabel('Conversión X_A [%]', 'FontWeight', 'bold');
legend('Lazo Abierto', 'Lazo Cerrado (PID)', 'Objetivo Nominal (73.0%)', 'Location', 'best');
grid on; xlim([0 60]);

% Subgráfica 3: Acción de control (Temperatura de chaqueta)
subplot(3, 1, 3);
plot(t_min, Tc_cl, '-m', 'LineWidth', 1.6); hold on;
plot(t_min, ones(size(t_min))*Tc_nom, ':k', 'LineWidth', 1.2);
xlabel('Tiempo de operación [min]', 'FontWeight', 'bold');
ylabel('Temp. Chaqueta T_c [K]', 'FontWeight', 'bold');
legend('Señal Manipulada T_c (PID)', 'T_c Nominal', 'Location', 'best');
grid on; xlim([0 60]);

% Exportar figura en alta resolución para README.md y Notion
print(fig, 'reactor_pid_resultados_es.png', '-dpng', '-r300');
