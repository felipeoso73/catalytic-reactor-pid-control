%% =========================================================================
%  CONTROLE PID DE TEMPERATURA EM REATOR CATALÍTICO BIFÁSICO (CSTR)
%  Processo: Valorização catalítica de resíduos agroindustriais a P = 25 bar
%  Conversão nominal: 73.0% a T_sp = 415.0 K
%  Instrumentação primária: Sensor RTD Pt100 Classe A (4 fios, poço AISI 316L)
%  Compatibilidade: MATLAB e GNU Octave
%  =========================================================================
clear; clc; close all;

%% 1. PARÂMETROS FÍSICO-QUÍMICOS E CONDIÇÕES NOMINAIS DE OPERAÇÃO
P_op     = 25.0;        % Pressão constante do sistema bifásico [bar]
V        = 2.5;         % Volume efetivo da suspensão líquido-catalisador [m^3]
F        = 0.005;       % Vazão volumétrica de alimentação [m^3/s]
tau      = V / F;       % Tempo espacial nominal [s] (500 s)

rho      = 850.0;       % Densidade média da mistura reacional [kg/m^3]
Cp       = 2800.0;      % Calor específico médio [J/(kg*K)]
dHr      = -95000.0;    % Entalpia de reação exotérmica [J/mol]
UA       = 1800.0;      % Coeficiente global de transferência de calor [W/K]

Ea       = 65000.0;     % Energia de ativação aparente (Ni-Mo/Al2O3) [J/mol]
R        = 8.314;       % Constante universal dos gases [J/(mol*K)]

T_sp     = 415.0;       % Temperatura de referência (Setpoint) [K]
X_nom    = 0.73;        % Eficiência de conversão nominal alvo [-] (73.0%)

% Calibração exata do fator pré-exponencial k0 para X_nom = 73.0% a 415 K
k_nom    = (X_nom / (1.0 - X_nom)) / tau;
k0       = k_nom * exp(Ea / (R * T_sp));

% Condições nominais de alimentação e regime estacionário
CA0_nom  = 1200.0;      % Concentração nominal do substrato residual [mol/m^3]
T0_nom   = 385.0;       % Temperatura nominal de alimentação [K]
CA_nom   = CA0_nom * (1.0 - X_nom);

% Balanço de energia estacionário para obter Tc nominal da camisa térmica
Q_conv   = rho * Cp * F * (T0_nom - T_sp);
Q_rxn    = (-dHr) * V * k_nom * CA_nom;
Tc_nom   = T_sp - (Q_conv + Q_rxn) / UA;

%% 2. INSTRUMENTAÇÃO PRIMÁRIA E PARÂMETROS DO CONTROLADOR PID
% Sensor industrial RTD Pt100 Classe A em poço termométrico AISI 316L
tau_sens = 8.0;         % Constante de tempo térmica do conjunto sensor [s]

% Sintonia do controlador PID (ação reversa sobre Tc)
Kc       = 4.5;         % Ganho proporcional [K/K]
tau_I    = 140.0;       % Tempo integral [s]
tau_D    = 12.0;        % Tempo derivativo [s]

% Limites físicos do fluido térmico da camisa (saturação)
Tc_min   = 320.0;       % Temperatura mínima de resfriamento [K]
Tc_max   = 430.0;       % Temperatura máxima de aquecimento [K]

%% 3. CONFIGURAÇÃO DA SIMULAÇÃO DINÂMICA E PERTURBAÇÕES DE CARGA
dt       = 1.0;         % Passo de amostragem e integração [s]
t_end    = 3600.0;      % Horizonte de simulação [s] (60 minutos)
t        = 0:dt:t_end;
N        = length(t);

% Vetores de armazenamento: Malha Aberta (OL) vs. Malha Fechada PID (CL)
T_ol     = zeros(1, N);
X_ol     = zeros(1, N);

T_cl     = zeros(1, N);
Tm_cl    = zeros(1, N);
X_cl     = zeros(1, N);
Tc_cl    = zeros(1, N);

% Estados iniciais em regime estacionário: [CA; T; Tm]
state_ol = [CA_nom; T_sp; T_sp];
state_cl = [CA_nom; T_sp; T_sp];

integral_err = 0.0;
prev_Tm      = T_sp;

% Equações diferenciais acopladas (balanço de massa, energia e atraso do sensor)
ode_rhs = @(s, CA0, T0, Tc) [ ...
    (F / V) * (CA0 - s(1)) - (k0 * exp(-Ea / (R * s(2)))) * s(1); ...
    (F / V) * (T0 - s(2)) + ((-dHr) / (rho * Cp)) * (k0 * exp(-Ea / (R * s(2)))) * s(1) ...
        - (UA / (rho * Cp * V)) * (s(2) - Tc); ...
    (s(2) - s(3)) / tau_sens ...
];

%% 4. LOOP DE INTEGRAÇÃO NUMÉRICA (RK4) E CONTROLE PID DISCRETO
for i = 1:N
    ti = t(i);
    
    % Perturbação 1 (t = 600 s): aumento degrau de +20% na concentração de resíduo
    if ti >= 600
        CA0 = 1.20 * CA0_nom;
    else
        CA0 = CA0_nom;
    end
    
    % Perturbação 2 (t = 1800 s): queda degrau de -10 K na temperatura de entrada
    if ti >= 1800
        T0 = T0_nom - 10.0;
    else
        T0 = T0_nom;
    end
    
    % --- Simulação em Malha Aberta (Tc fixa em Tc_nom) ---
    T_ol(i) = state_ol(2);
    X_ol(i) = ((CA0 - state_ol(1)) / CA0) * 100.0;
    
    k1 = ode_rhs(state_ol, CA0, T0, Tc_nom);
    k2 = ode_rhs(state_ol + 0.5*dt*k1, CA0, T0, Tc_nom);
    k3 = ode_rhs(state_ol + 0.5*dt*k2, CA0, T0, Tc_nom);
    k4 = ode_rhs(state_ol + dt*k3, CA0, T0, Tc_nom);
    state_ol = state_ol + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
    
    % --- Simulação em Malha Fechada (Controle PID com RTD Pt100) ---
    T_cl(i)  = state_cl(2);
    Tm_cl(i) = state_cl(3);
    X_cl(i)  = ((CA0 - state_cl(1)) / CA0) * 100.0;
    
    err      = T_sp - Tm_cl(i);
    dTm_dt   = (Tm_cl(i) - prev_Tm) / dt;
    prev_Tm  = Tm_cl(i);
    
    % Lei de controle PID (derivada na variável medida para evitar kick derivativo)
    u_unsat  = Tc_nom + Kc * (err + (1.0 / tau_I) * integral_err - tau_D * dTm_dt);
    Tc_cmd   = min(max(u_unsat, Tc_min), Tc_max);
    
    % Integração condicional (Anti-Windup)
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

%% 5. CÁLCULO DE INDICADORES DE DESEMPENHO E RESUMO NO CONSOLE
IAE_ol = trapz(t, abs(T_sp - T_ol));
IAE_cl = trapz(t, abs(T_sp - T_cl));

fprintf('=================================================================\n');
fprintf(' RESUMO DE DESEMPENHO - REATOR CATALÍTICO BIFÁSICO (P = %.1f bar)\n', P_op);
fprintf('=================================================================\n');
fprintf(' Conversão nominal inicial        : %.2f %%\n', X_cl(1));
fprintf(' Temperatura máxima (Malha Aberta): %.2f K\n', max(T_ol));
fprintf(' Temperatura máxima (Controle PID): %.2f K\n', max(T_cl));
fprintf(' Conversão final (Malha Aberta)   : %.2f %%\n', X_ol(end));
fprintf(' Conversão final (Controle PID)   : %.2f %%\n', X_cl(end));
fprintf(' Erro Integral Absoluto IAE (OL)  : %.1f K*s\n', IAE_ol);
fprintf(' Erro Integral Absoluto IAE (PID) : %.1f K*s (Redução: %.1f%%)\n', ...
        IAE_cl, (1 - IAE_cl/IAE_ol)*100);
fprintf('=================================================================\n');

%% 6. GERAÇÃO DE GRÁFICOS COMPARATIVOS E EXPORTAÇÃO PARA PORTFÓLIO
t_min = t / 60.0;
fig = figure('Color', 'w', 'Position', [100, 100, 950, 720]);

% Subgráfico 1: Temperatura do reator
subplot(3, 1, 1);
plot(t_min, T_ol, '--r', 'LineWidth', 1.6); hold on;
plot(t_min, T_cl, '-b',  'LineWidth', 1.8);
plot(t_min, ones(size(t_min))*T_sp, ':k', 'LineWidth', 1.2);
ylabel('Temperatura [K]', 'FontWeight', 'bold');
title('Controle PID de Temperatura em Reator Catalítico Bifásico (P = 25 bar)', 'FontSize', 11);
legend('Malha Aberta', 'Malha Fechada (PID + Pt100)', 'Setpoint (415 K)', 'Location', 'best');
grid on; xlim([0 60]);

% Subgráfico 2: Eficiência de conversão
subplot(3, 1, 2);
plot(t_min, X_ol, '--r', 'LineWidth', 1.6); hold on;
plot(t_min, X_cl, '-g',  'LineWidth', 1.8);
plot(t_min, ones(size(t_min))*73.0, ':k', 'LineWidth', 1.2);
ylabel('Conversão X_A [%]', 'FontWeight', 'bold');
legend('Malha Aberta', 'Malha Fechada (PID)', 'Alvo Nominal (73.0%)', 'Location', 'best');
grid on; xlim([0 60]);

% Subgráfico 3: Ação de controle (Temperatura da camisa)
subplot(3, 1, 3);
plot(t_min, Tc_cl, '-m', 'LineWidth', 1.6); hold on;
plot(t_min, ones(size(t_min))*Tc_nom, ':k', 'LineWidth', 1.2);
xlabel('Tempo de operação [min]', 'FontWeight', 'bold');
ylabel('Temp. Camisa T_c [K]', 'FontWeight', 'bold');
legend('Sinal Manipulado T_c (PID)', 'T_c Nominal', 'Location', 'best');
grid on; xlim([0 60]);

% Exportar figura em alta resolução para README.md e Notion
print(fig, 'reactor_pid_resultados_pt.png', '-dpng', '-r300');
