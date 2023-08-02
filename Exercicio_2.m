%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equipe:                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ademar A. Santos Jr.    %
% Leonardo Pessôa         %
% Hebert Crispim          %
% Edgley Carvalho         %
% Marcus Vinícius Pereira %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Descrição da Atividade 
% Exercício 2:
% Trata-se de um circuito magnético com dois enrolamentos que estão
% acoplados magneticamente. São fornecidos dados acerca da geometria do
% problema, dos terminais ponto, do número de espiras dos enrolamentos, da
% frequência da rede, bem como da curva de saturação do material
% ferromagnético do núcleo.


% São solicitadas 4 simulações, apresentando correntes normalizadas:
% 1) Carga: Resistor 2 Ohms,    Rede: tensão senoidal 60 Hz e 30 Vrms;
% 2) Carga: Resistor 2000 Ohms, Rede: tensão senoidal 60 Hz e 30 Vrms; 
% 3) Carga: Resistor 2 Ohms,    Rede: tensão senoidal 60 Hz e 90 Vrms;
% 4) Carga: Resistor 2000 Ohms, Rede: tensão senoidal 60 Hz e 90 Vrms; 


% Além disso, são solicitados textos discussivos sobre o motivo das
% senóides apresentarem distorção ou não e sobre a ausência de perdas dos
% cenários simulados.

%% Configurações do Arquivo
clear; close all; clc;      % Limpeza do terminal, das figuras e das variáveis
format longE;               % Variáveis com 15 casas decimais e em notação científica

% Abrir arquivo auxiliar
Tabela_BxH_Nucleo;          % Arquivo contendo dados da curva BxH do ímã

%% Dados da questão (SI)
% Dados de geometria do sistema
Ac = 2e-4;                  % Área do Núcleo
lc = 30e-2;                 % Comprimento médio total do núcleo

% Dados das bobinas
N1 = 1000;                  % Número de espiras do enrolamento 1
N2 = 2000;                  % Número de espiras do enrolamento 2

% Dados da rede elétrica
f = 60;                     % Frequência da rede elétrica

%% Parâmetros da simulação
h = 1e-4;                   % Passo de cálculo
t0 = 0;                     % Tempo inicial
tf = 1;                     % Tempo final
n = (tf - t0)/h;            % Número de pontos
k = 0;                      % Posição dos vetores de salvamento

%% Parâmetros Adicionais
w = 2*pi*f;                 % Frequência angular da rede

% Densidade de fluxo do núcleo
Bc = linspace(0,1.8,n);     % Valores arbitrados de densidade de fluxo no núcleo

% Intensidade de campo mangético do núcleo
Hc = Hfun_real(Bc);         % Intensidade de campo no núcleo

% Permeabilidade e relutância do núcleo para o caso ideal
mi_ideal = 1e12;            % Permeabilidade ideal do ferro (~inf)
Rl_ideal = lc/(mi_ideal*Ac);% Relutância ideal do ferro (~0)

%% Condições do Sistema
% Tensão eficaz da rede
Vef = 30;                   % Valor eficaz da tensão da rede elétrica

% Resistência da carga
R = 2;                      % Resistência da carga

%% Valores máximos
% Cálculo da tensão máxima
Vf_max = Vef*sqrt(2);       % Tensão máxima da rede elétrica
Ve_max = (N2/N1)*Vf_max;    % Tensão máxima na carga

% Cálculo do fluxo máximo
fluxo_max = Vf_max/(N1*w);  % Fluxo máximo

% Correntes máximas para o caso ideal
i1_ideal_max = (fluxo_max*Rl_ideal)/N1 + (Vf_max/R)*(N2/N1)^2; % Corrente ideal máxima da bobina 1
i2_ideal_max = (N1/N2)*i1_ideal_max; % Corrente máxima da bobina 2

%% Condições Iniciais
lambda1 = 0;                % Fluxo concatenado da bobina 1
lambda20 = 0;               % Fluxo concatenado da bobina 2 na posição n-1 (Para derivação)

%% Início do Looping
for t = t0:h:tf 
    %% Simulação
    % Tensão da rede elétrica
    Vf = Vf_max*cos(w*t);
        
    % Fluxos
    lambda1 = lambda1 + h*Vf;       % Fluxo concatenado - Bobina 1
    fluxo = lambda1/N1;             % Fluxo
    lambda2 = fluxo*N2;             % Fluxo concatenado - Bobina 2
    
    % Tensão sobre a carga
    Ve = (lambda2-lambda20)/h;      
    
    lambda20 = lambda2;             % Avanço do fluxo concatenado na posição n-1
    
    % Características instantâneas do núcleos
    B_inst = fluxo/(Ac);            % B devido à existência do λ
    H_inst = Hfun_real(B_inst);     % H real devido B
    mi_inst = B_inst/H_inst;        % Permeabilidade instantânea
    Rl_inst = lc/(mi_inst*Ac);      % Relutância instantânea do núcleo
    
    % Correntes
    i1_real  = (fluxo*Rl_inst)/N1 + (Vf/R)*(N2/N1)^2;
    i1_ideal = (fluxo*Rl_ideal)/N1 + (Vf/R)*(N2/N1)^2;
    
    i2_real  = Ve/R;
    i2_ideal = (N1/N2)*i1_ideal;
    
    % Correntes máximas
    i1_real_max  = (fluxo_max*Rl_inst)/N1 + (Vf_max/R)*(N2/N1)^2;  % Corrente real máxima da bobina 1
    i2_real_max  = Ve_max/R;
    
    %% Salvamento das variáveis
    
    % Progressão do índice de salvamento
    k = k+1;
    
    % Tempo da simulação
    Ts(k) = t;
    
    % Tensões
    Vfs(k) = Vf;
    Ves(k) = Ve;
    
    % Fluxos
    fluxos(k) = fluxo;
    lambda1s(k) = lambda1;
    lambda2s(k) = lambda2;
    
    % Correntes
    % i1
    i1_reals(k) = i1_real;
    i1_ideals(k) = i1_ideal;
    % i2
    i2_reals(k) = i2_real;
    i2_ideals(k) = i2_ideal;
    
    % Corrente máxima
    i1_real_maxs(k) = i1_real_max;
    i2_real_maxs(k) = i2_real_max;
    
end
% Correntes Normalizadas
% i1
i1N_reals = i1_reals./i1_real_maxs;
i1N_ideals = i1_ideals./i1_ideal_max;
% i2
i2N_reals = i2_reals./i2_real_maxs;
i2N_ideals = i2_ideals./i2_ideal_max;

%% Plots
% CURVA BxH
% figure('Name', 'Curva BxH do núcleo')
% plot(Hc,Bc,'k-','LineWidth',1),zoom
% title('Curva BxH')
% xlabel('H [ A/m ]','fontweight','bold','Fontsize',12)
% ylabel('B [ Wb/m² ]','fontweight','bold','Fontsize',12)
% axis([0 1500 -0.18 2])
% grid minor
% set(gca, 'XScale', 'log') % Escala log

% % CURVA λxI1
% I1 = Hc.*lc/N1; % Corrente de Magnetização 
% Lb = Bc.*Ac;
% figure('Name', 'Curva λxi: 1')
% plot(I1,Lb,'k-','LineWidth',1),zoom
% title('Curva λxi: 1')
% xlabel('$i_1$ [ A$e$ ]','fontweight','bold','Fontsize',14,'interpreter','latex')
% ylabel('$\lambda_1$ [ Wb$e$ ]','fontweight','bold','Fontsize',14,'interpreter','latex')
% axis([0 1.5 0 5e-4])
% grid minor

% % TENSÕES
% figure('Name','Tensão de entrada e saída')
% plot(Ts,Vfs,'r-',Ts,Ves,'b-'),zoom
% title('Tensão de entrada (Rede) e saída (Carga)')
% legend('$v_{Rede}$','$v_{Carga}$','Fontsize',16,'interpreter','latex')
% ylabel('Tens\~{a}o [ V ]','fontweight','bold','Fontsize',14,'interpreter','latex')
% xlabel('Tempo [ $s$ ]','fontweight','bold','Fontsize',14,'interpreter','latex')
% grid minor

% % FLUXO
% figure('Name','Fluxo')
% plot(Ts,fluxos),zoom
% title('Fluxo no ferro')
% ylabel('Fluxo [ $Wb$ ]','fontweight','bold','Fontsize',14,'interpreter','latex')
% xlabel('Tempo [ $s$ ]','fontweight','bold','Fontsize',14,'interpreter','latex')
% axis([3*tf/4 tf -5e-4 5e-4])
% grid minor

% CORRENTES NÃO NORMALIZADAS
figure('Name','Corrente 1 e 2')
subplot(2,1,1)
plot(Ts,i1_reals,'r-',Ts,i1_ideals,'k--'),zoom
title('Corrente i1')
legend('$i_1$ Real','$i_1$ Ideal','Fontsize',16,'interpreter','latex')
ylabel('$\bf{I [ A ]}$','Fontsize',14,'interpreter','latex')
xlabel('$\bf{t [ s ]}$','Fontsize',14,'interpreter','latex')
axis([3*tf/4 tf -0.6 0.6])
grid minor
subplot(2,1,2)
plot(Ts,i2_reals,'b-',Ts,i2_ideals,'k--'),zoom
title('Corrente i2')
legend('$i_2$ Real','$i_2$ Ideal','Fontsize',16,'interpreter','latex')
ylabel('$\bf{I [ A ]}$','Fontsize',14,'interpreter','latex')
xlabel('$\bf{t [ s ]}$','Fontsize',14,'interpreter','latex')
grid minor
axis([3*tf/4 tf -0.6 0.6])

% CORRENTES NORMALIZADAS
figure('Name','Corrente 1 e 2 normalizadas')
subplot(2,1,1)
plot(Ts,i1N_reals,'r-',Ts,i1N_ideals,'k--'),zoom
title('Corrente i1 Normalizada')
legend('$i_1$ Real','$i_1$ Ideal','Fontsize',16,'interpreter','latex')
ylabel('$\bf{I [ A ]}$','Fontsize',14,'interpreter','latex')
xlabel('$\bf{t [ s ]}$','Fontsize',14,'interpreter','latex')
axis([3*tf/4 tf -1.1 1.1])
grid minor
subplot(2,1,2)
plot(Ts,i2N_reals,'b-',Ts,i2N_ideals,'k--'),zoom
title('Corrente i2 Normalizada')
legend('$i_2$ Real','$i_2$ Ideal','Fontsize',16,'interpreter','latex')
ylabel('$\bf{I [ A ]}$','Fontsize',14,'interpreter','latex')
xlabel('$\bf{t [ s ]}$','Fontsize',14,'interpreter','latex')
grid minor
axis([3*tf/4 tf -1.1 1.1])