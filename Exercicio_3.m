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
% Exercício 3:
% Consiste de um circuito magnético contendo uma peça móvel e
% um íma permanente. É fornecido um trecho do laço de Histerese do ímã.


% São consideradas duas modelagens para o laço de histerese: 
% 1) Realização de um fitting da curva BxH fornecida;
% 2) Apoximação Linear da curva BxH fornecida.


% São solicitadas 3 atividades:

% a) Apresentação do gráfico da força magnética que surge na peça móvel em
% função do seu deslocamento;

% b) Instalação de uma bobina de 1500e:
% Parte 1: Apresentação do gráfico da densidade de fluxo no entreferro em
% função da corrente aplicada na bobina;
% Parte 2: Apresentação da força magnética na peça móvel em função do seu
% deslocamento no cenário em que a bobina está sendo excitada por uma
% corrente que coloque o sistema no ponto de operação normal do gráfico no
% deslocamento intermediário.

%% Configurações do Arquivo
clear; close all; clc;      % Limpeza do terminal, das figuras e das variáveis
format longE;               % Variáveis com 15 casas decimais e em notação científica

% Abrir arquivo auxiliar
Tabela_BxH_Ima;             % Arquivo contendo dados da curva BxH do ímã 

%% Dados da questão (SI)
% Dados da bobina fictícia (Questão a)
Nf = 1000;                  % Número de espiras fictício
If = 1;                     % Corrente fictícia

% Dados da bobina inserida (Questão b)
N = 1.5e3;                  % Número de espiras (Questão 2)

% Dados de geometria do sistema
g = 2e-3;                   % Comprimento do entreferro
Wm = 2e-2;                  % Largura do ímã
Wg = 2.5e-2;                % Largura do entreferro
d = 1e-2;                   % Comprimento do ímã
D = 3e-2;                   % Largura do núcleo

Ai = Wm*D;                  % Área do ímã

%% Parâmetros da simulação
h = 1e-2;                   % Passo de cálculo
t0 = 0;                     % Tempo inicial
tf = 1;                     % Tempo final
n = (tf - t0)/h;            % Número de pontos
k = 0;                      % Posição dos vetores de salvamento

%% Parâmetros Adicionais
% Dados da geometria móvel
x  = linspace(0,Wg-1e-10,n);% Deslocamento do êmbolo
Ag = D.*(Wg-x);             % Área do êmbolo disponível para passagem do fluxo  

% Permeabilidade do ar (Entreferro)
mi0 = 4*pi*10^(-7);         % Permeabilidade magnética no ar

% Permeabilidade do núcleo
mi = inf;                   % Permeabilidade magnética no núcleo

% Densidades de fluxo
Bi = linspace(0,1.2,n);     % Valores arbitrados para a densidade de fluxo no ímã
Bg = (2*Wm/Wg).*Bi;         % Densidade de fluxo no entreferro

% Intensidade de campo mangético do ímã (Duas modelagens)
Hireal = Hfun_real(Bi);     % Intensidade de campo do ímã - Caso real
Hilinear = Hfun_linear(Bi); % Intensidade de campo do ímã - Caso linear

% Corrente de coercitividade
Ifcoerc = (Hcoerc*d)/Nf;    % Corrente fictícia que anula o fluxo do ímã

% Fluxo concatenado
lambdaf = Nf*Ai.*Bi;        % Fluxo concatenado na bobina fictícia Nf
lambda = N*Ai*Bi;           % Fluxo concatenado na bobina (Questão b)

%% Questão 1
% Looping de simulação progredindo o deslocamento da peça móvel
for dx = x
    % Cálculo da corrente fictícia em função do fluxo concatenado fictício
    % Pela lei de Ampère:
    if_real     = (d/Nf).*Hireal   + (2*g.*lambdaf)./(D*(Nf^2)*mi0*(Wg-dx)); % Curva BxH real
    if_linear   = (d/Nf).*Hilinear + (2*g.*lambdaf)./(D*(Nf^2)*mi0*(Wg-dx)); % Aproximação Linear
    
    % Definição da função Corrente fictícia vs Fluxo Concatenado fictício
    fifun_real   = @(I0) spline(if_real  , lambdaf, I0); % Curva BxH real
    fifun_linear = @(I0) spline(if_linear, lambdaf, I0); % Aproximação Linear
    
    % Cálculo da coenergia do sistema (Área abaixo da curva λf x If)
    cW_real   = integral(fifun_real,  Ifcoerc,0); % Curva BxH real
    cW_linear = integral(fifun_linear,Ifcoerc,0); % Aproximação Linear
    
    %% Salvamento das Variáveis
    % progressão do índice do vetor de salvamento
    k = k+1;

    % Deslocamento
    xs(k)  = dx; 
    
    % Coenergia
    cW_reals(k) = cW_real;     % Curva BxH real
    cW_linears(k) = cW_linear; % Aproximação Linear
end
% Fim do Looping

% Cálculo da força magnética (Derivada da coenergia no deslocamento)
F_reals   = diff(cW_reals)./diff(xs); % Curva BxH real
F_linears = diff(cW_linears)./diff(xs); % Aproximação Linear

% Ajuste do array de posição (O processo de diferenciação gera um array com n-1 posições)
xfs = xs;                  % Array de deslocamento relacionado à força
xfs(1) = [];               % Exclusão do primeiro elemento

%% Plots - Questão a

%  Coenergia vs Deslocamento
% figure('Name', 'Coenergia vs Deslocamento do êmbolo - a')
% plot(xs.*100,cW_reals,'r-',xs.*100,cW_linears,'b-','LineWidth',2),zoom
% title('Coenergia vs Deslocamento - Questão a)','Fontsize',16)
% xlabel('\bf{X [ cm ]}','Fontsize',12,'interpreter','latex')
% ylabel("\bf{W' [ J ]}",'Fontsize',12,'interpreter','latex')
% legend("$W_{Real}'$","$W_{Linear}'$",'Fontsize',16,'interpreter','latex')
% grid minor
% Força vs Deslocamento 
figure('Name', 'Força vs Deslocamento do êmbolo - Questão a)')
plot(xfs.*100,F_reals,'r-',xfs.*100,F_linears,'b-','LineWidth',2),zoom
title('Força vs Deslocamento - Questão a)','Fontsize',16)
xlabel('\bf{X [ cm ]}','Fontsize',12,'interpreter','latex')
ylabel('\bf{F [ N ]}','Fontsize',12,'interpreter','latex')
legend('$F_{Real}$','$F_{Linear}$','Fontsize',16,'interpreter','latex')
grid minor

%% Questão b
% Parte 1

% Cálculo da corrente na bobina para a peça móvel na posição intermediária (x = wg/2)
% Pela lei de Ampère:
I_real   = Hireal.*(d/N)   + Bg.*(2*g/(mi0*N));
I_linear = Hilinear.*(d/N) + Bg.*(2*g/(mi0*N));

% Parte 2
% Corrente no ponto de operação para a peça móvel na posição intermediária (x = wg/2)
% Pela lei de Ampère:
Icte_real =   Hopn_real*(d/N)     + Bopn*(2*Wm/Wg)*(2*g/(mi0*N));
Icte_linear = Hopn_linear*(d/N)   + Bopn*(2*Wm/Wg)*(2*g/(mi0*N));

% Reinicialização da variável de salvamento
k = 0; 

% Looping de simulação progredindo o deslocamento da peça móvel
for dx = x
    % Cálculo da corrente em função do fluxo concatenado
    % Pela lei de Ampère:
    i_real     = (d/N).*Hireal   + (2*g.*lambda)./((N^2)*mi0*D*(Wg-dx));
    i_linear   = (d/N).*Hilinear + (2*g.*lambda)./((N^2)*mi0*D*(Wg-dx));
    
    % Definição da função Corrente vs Fluxo Concatenado
    fifun_real1   = @(I0) spline(i_real  , lambda, I0);
    fifun_linear1 = @(I0) spline(i_linear, lambda, I0);
    
    % Cálculo da coenergia do sistema (Área abaixo da curva λ x I)
    cW_real1   = integral(fifun_real1,  Ifcoerc,Icte_real);
    cW_linear1 = integral(fifun_linear1,Ifcoerc,Icte_linear);
    
    %% Salvamento das Variáveis
    k = k+1;

    % Deslocamento
    x1s(k)  = dx; 
    
    % Coenergia
    cW_real1s(k) = cW_real1;
    cW_linear1s(k) = cW_linear1;
    
end
% Fim do Looping

% Cálculo da força magnética (Derivada da coenergia no deslocamento)
F_real1s   = diff(cW_real1s)./diff(x1s); % Curva BxH real
F_linear1s = diff(cW_linear1s)./diff(x1s); % Aproximação Linear

xf1s = x1s;                 % Array de deslocamento relacionado à força
xf1s(1) = [];               % Exclusão do primeiro elemento


%% Plots - Questão b
% Parte 1
% Densidade de fluxo no entreferro vs Corrente de excitação
figure('Name', 'Densidade de fluxo no entreferro vs Deslocamento do êmbolo - Questão b1)')
plot(I_real,Bg,'r-',I_linear,Bg,'b-','LineWidth',2),zoom
title('Densidade de fluxo no entreferro vs Deslocamento - Questão b1)','Fontsize',16)
xlabel('\bf{I [ A ]}','Fontsize',12,'interpreter','latex')
ylabel('$\bf{B_{g} [ \frac{Wb}{m^2} ]}$','Fontsize',12,'interpreter','latex')
legend('$I_{Real}$','$I_{Linear}$','Fontsize',16,'interpreter','latex')
grid minor

% Parte 2
% Coenergia vs Deslocamento 
% figure('Name', 'Coenergia vs Deslocamento do êmbolo - Questão b2)')
% plot(x1s.*100,cW_real1s,'r-',x1s.*100,cW_linear1s,'b-','LineWidth',2),zoom
% title('Coenergia vs Deslocamento - Questão b2)','Fontsize',16)
% xlabel('\bf{X [ cm ]}','Fontsize',12,'interpreter','latex')
% ylabel("\bf{W' [ J ]}",'Fontsize',12,'interpreter','latex')
% legend("$W_{Real}'$","$W_{Linear}'$",'Fontsize',16,'interpreter','latex')
% grid minor
% Força vs Deslocamento
figure('Name', 'Força vs Deslocamento do êmbolo - Questão b2)')
plot(xf1s.*100,F_real1s,'r-',xf1s.*100,F_linear1s,'b-','LineWidth',2),zoom
title('Força vs Deslocamento  - Questão b2)','Fontsize',16)
xlabel('\bf{X [ cm ]}','fontweight','bold','Fontsize',12,'interpreter','latex')
ylabel('\bf{F [ N ]}','fontweight','bold','Fontsize',12,'interpreter','latex')
legend('$F_{Real}$','$F_{Linear}$','Fontsize',16,'interpreter','latex')
grid minor