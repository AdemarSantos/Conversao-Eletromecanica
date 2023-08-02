%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equipe:                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ademar A. Santos Jr.    %
% Leonardo Pessôa         %
% Hebert Crispim          %
% Edgley Carvalho         %
% Marcus Vinícius Pereira %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Configurações do Arquivo

clear; close all; clc; % Limpeza do terminal, das figuras e variáveis
format longE; % variáveis com 15 casas decimais e em notação científica

% Abrir arquivo auxiliar
Tabela_BxH_Nucleo; % Arquivo contendo dados da curva BxH do ímã

%% Descrição da Atividade
% Exercício 1 - Questão 4 
% Apresente um gráfico do tempo necessário para o elemento móvel sair da posição x=d/2 até a
% posição x=0 em função da corrente aplicada na bobina. Seu gráfico deve compreender correntes
% de zero até o valor limite calculado no item anterior.

% ***** TEMPO DE SIMULAÇÃO DE CERCA DE 90 SEGUNDOS *****

%% Dados da questão
N = 1e3;                    % Número de espiras
g = 2e-3;                   % Espaçamento do entreferro
d = 4e-2;                   % Largura do núcleo (CTE)
l = 4e-2;                   % Espessura do núcleo (CTE)
lc = 70e-2;                 % Comprimento total pelo caminho médio do núcleo
mi0 = 4*pi*10^(-7);         % Permeabilidade magnética no vácuo (apx ar)

%% Parâmetros da simulação
h = 1e-2;                   % Passo de cálculo
t0 = 0;                     % Tempo inicial
tf = 1;                     % Tempo final
n = (tf - t0)/h;            % Número de pontos
k = 0;                      % Posição dos vetores de salvamento

%% Parâmetros Adicionais
% Tomando como corrente de confiabilidade o último valor de corrente que pode ser obtida através da tabela BxH,
% temos que IF = 19.73 A, ou seja, teremos valores confiáveis para um valor de I < IF. Dessa forma, foi utilizado 
% IF  = 19.7 A como corrente máxima para realização das atividades.

If = 19.7;                  % Corrente de confiabilidade (valor aproximado)
I  = linspace(1e-6,If,n);   % Array de corrente para integração
x  = linspace(0,d-1e-10,n); % Deslocamento do êmbolo
Ac = l*d;                   % Área do núcleo
Ag = l*(d-x);               % Área do êmbolo disponível para passagem do fluxo
temporealS = zeros(n,0);    % Variável para alocação do vetor tempo real
tempoidealS = zeros(n,0);   % Variável para alocação do vetor tempo ideal
tempolinearS = zeros(n,0);  % Variável para alocação do vetor tempo linear

Bc = linspace(0,1.8,n);     % Valores arbitrados de densidade de fluxos

Hc_real = Hfun_real(Bc);        % Intensidade de campo no núcleo - Caso real
Hc_linear = Hfun_linear(Bc);    % Intensidade de campo no núcleo - Caso linear
Hc_ideal = 0.*Bc;           % Intensidade de campo no núcleo - Caso ideal

% Fluxos concatenados para todos os valores de Bc
lambda = N*Ac.*Bc;

%% Iniciando Loop 1
for i = 1:1:n % Loop geral que calcula o valor do tempo pra uma dada corrente e salva este valor para as nossas n correntes
     for dx = x % Loop de inicialização e cálculo de nossas funções e variáveis
         
        % definição de nossa função corrente
        i_real   = (Hc_real*lc)/N + 2*g.*lambda./(N^2*l*mi0*(d-dx));
        i_linear = (Hc_linear*lc)/N + 2*g.*lambda./(N^2*l*mi0*(d-dx));
        i_ideal  = (Hc_ideal*lc)/N + 2*g.*lambda./(N^2*l*mi0*(d-dx));

        % Função Corrente vs Fluxo concatenado'
        fifun_real   = @(I0) spline(i_real, lambda, I0);
        fifun_linear = @(I0) spline(i_linear, lambda, I0);
        fifun_ideal  = @(I0) spline(i_ideal, lambda, I0);

        % Coenergia
        cW_real   = integral(fifun_real,0,I(i));    
        cW_linear = integral(fifun_linear,0,I(i));
        cW_ideal  = integral(fifun_ideal,0,I(i));

        %% Salvamento
        k = k+1;

        for m = 1:1:100
            i_reals(m,k) = i_real(m);
        end

        % Deslocamento
        xs(k)  = dx; 

        % Coenergia
        cW_reals(k) = cW_real;
        cW_linears(k) = cW_linear;
        cW_ideals(k) = cW_ideal;

     end
    
    % Calculo da força
    F_reals   = diff(cW_reals)./diff(xs);
    F_linears = diff(cW_linears)./diff(xs);
    F_ideals  = diff(cW_ideals)./diff(xs);

    %% Inicialização das Variáveis
    cW0_real = 0;
    cW0_linear = 0;
    cW0_ideal = 0;

    xfs = xs;
    xfs(1) = []; % Apagando o primeiro elemento para evitar problemas de convergência
    
    % Assumindo que a massa da peça móvel seja de 1 kg :
    m = 1; 
    % Logo, como F = m*a, nossa aceleração será igual a força
    
    % Aceleração em função da posição
    a_reals   = @(X0) spline(xfs, F_reals./m, X0);
    a_linears   = @(X0) spline(xfs, F_linears./m, X0);
    a_ideals   = @(X0) spline(xfs, F_ideals./m, X0);

    x = linspace(d/2,1e-12, n);
    k = 0;

    % Calculando a velocidade da peça
    for dx = x
        
        % Integral da função aceleração no delocamento de d/2 até x
        wreal = integral(a_reals,d/2,dx);
        wideal = integral(a_ideals,d/2,dx); 
        wlinear = integral(a_linears,d/2,dx); 

        k = k+1;
        % Salvamento dos valores de velocidade
        wreals(k) = wreal;
        wideals(k) = wideal; 
        wlinears(k) = wlinear; 
    end
    
    x(1) = [];   % Eliminando o primeiro ponto para evitar problemas de convergência (1/v, v=!0)
    Velreal = 2.*sqrt(wreals); % Velocidade real
    Velideal = 2.*sqrt(wideals); % Velocidade ideal
    Vellinear = 2.*sqrt(wlinears); % Velocidade linear
    
    % Eliminando o primeiro ponto para evitar problemas de convergência (1/v, v=!0)
    Velreal(1) = []; 
    Vellinear(1) = []; 
    Velideal(1) = []; 
    
    Vel_reals = @(X0) spline(x, Velreal, X0); % Função v(x)
    invVel_reals = @(X0) spline(x,1./Vel_reals(x), X0); % Função 1/v(x)
    
    Vel_ideals = @(X0) spline(x, Velideal, X0); % Função v(x)
    invVel_ideals = @(X0) spline(x,1./Vel_ideals(x), X0); % Função 1/v(x)
    
    Vel_linears = @(X0) spline(x, Vellinear, X0); % Função v(x)
    invVel_linears = @(X0) spline(x,1./Vel_linears(x), X0); % Função 1/v(x)
    
    % Cálculo do tempo real
    temporeal = abs(integral(invVel_reals,x(1),x(length(x)))); 
    temporealS(i) = temporeal; % salvamento do tempo deste loop no nosso array
    
    % Cálculo do tempo ideal
    tempoideal = abs(integral(invVel_ideals,x(1),x(length(x)))); 
    tempoidealS(i) = tempoideal; % salvamento do tempo deste loop no nosso array
    
    % Cálculo do tempo linear
    tempolinear = abs(integral(invVel_linears,x(1),x(length(x)))); 
    tempolinearS(i) = tempolinear; % salvamento do tempo deste loop no nosso array
    
end 

% O tempo ficou negativo, acredito que fruto de alguma troca de sinal que
% ficou incorreta, por isso tirei o absoluto (A força é negativa, a acelera
% tbm vai ser, a integral da aceleração de d/2 até 0 faz um resultado
% positivo (velocidade), a integral da velocidade (positiva) de d/2 até
% zero, dá um valor negativo...

% Tempo Real
figure('Name', 'Tempo para se deslocar de d/2 a 0 - Real')
title('Tempo para se deslocar de d/2 a 0','Fontsize',16)
plot(I(2:1:n), temporealS(2:1:n),'LineWidth',3),zoom
legend('Tempo - Força Real')
ylabel('Tempo [ s ]','fontweight','bold','Fontsize',12)
xlabel('I [ A ]','fontweight','bold','Fontsize',12)
grid minor 

% Tempo Ideal
figure('Name', 'Tempo para se deslocar de d/2 a 0 - Ideal')
title('Tempo para se deslocar de d/2 a 0','Fontsize',16)
plot(I(2:1:n), tempoidealS(2:1:n),'LineWidth',3),zoom
legend('Tempo - Força Ideal')
ylabel('Tempo [ s ]','fontweight','bold','Fontsize',12)
xlabel('I [ A ]','fontweight','bold','Fontsize',12)
grid minor 

% Tempo Linear
figure('Name', 'Tempo para se deslocar de d/2 a 0 - Linear')
title('Tempo para se deslocar de d/2 a 0','Fontsize',16)
plot(I(2:1:n), tempolinearS(2:1:n),'LineWidth',3),zoom
legend('Tempo - Força Linear')
ylabel('Tempo [ s ]','fontweight','bold','Fontsize',12)
xlabel('I [ A ]','fontweight','bold','Fontsize',12)
grid minor 