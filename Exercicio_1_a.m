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
% Exercício 1:

% Simular força magnética que surge em um êmbolo para o caso
% ideal, linear e real do núcleo.

% 1) Considerando uma corrente aplicada na bobina de 1 A,
% apresente um gráfico para a força magnética que surge
% no sistema em função do comprimento do entreferro nas condições a seguir:

%   a)O núcleo como sendo ideal;
%   b)Uma aproximação linear da característica do núcleo;
%   c)A característica real do núcleo.

% 2) Repita o exercício anterior para uma corrente de 5 A e para uma corrente de 15 A.

% 3) Considerando os dados disponíveis, qual seria a maior corrente que
% poderia ser levar em consideração de modo que se tenha confiabilidade no cálculo da 
% força magnética em função da posição do elemento móvel.

% Para a questão quatro, utilizaremos outro arquivo

%% Configurações do Arquivo

clear; close all; clc; % Limpeza do terminal, das figuras e variáveis
format longE; % variáveis com 15 casas decimais e em notação científica

% Abrir arquivo auxiliar
Tabela_BxH_Nucleo; % Arquivo contendo dados da curva BxH do ímã


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

%% Parâmetros Adicionais
x  = linspace(0,d-1e-10,n); % Deslocamento do êmbolo
Ac = l*d;                   % Área do núcleo
Ag = l*(d-x);               % Área do êmbolo disponível para passagem do fluxo  


Bc = linspace(0,1.8,n);     % Valores arbitrados de densidade de fluxos

Hc_real = Hfun_real(Bc);        % Intensidade de campo no núcleo - Caso real
Hc_linear = Hfun_linear(Bc);    % Intensidade de campo no núcleo - Caso linear
Hc_ideal = 0.*Bc;           % Intensidade de campo no núcleo - Caso ideal

% Fluxos concatenados para todos os valores de Bc
lambda = N*Ac.*Bc;

%% Início do Looping
prompt = input('Insira um valor de corrente: '); % Solicita um valor de corrente para o usuário 
If = prompt;                % Corrente de excitação
I  = linspace(0,If,n);      % Array de corrente para integração
while(prompt)
    % Inicializa as variáveis do loop
    k = 0;                      % Posição dos vetores de salvamento
    cW0_real = 0;
    cW0_linear = 0;
    cW0_ideal = 0;
    F_reals   = 0;
    F_linears = 0;
    F_ideals  = 0;
    for dx = x
        i_real   = (Hc_real*lc)/N + 2*g.*lambda./(N^2*l*mi0*(d-dx));
        i_linear = (Hc_linear*lc)/N + 2*g.*lambda./(N^2*l*mi0*(d-dx));
        i_ideal  = (Hc_ideal*lc)/N + 2*g.*lambda./(N^2*l*mi0*(d-dx));

        % Função Corrente vs Fluxo concatenado'
        fifun_real   = @(I0) spline(i_real, lambda, I0);
        fifun_linear = @(I0) spline(i_linear, lambda, I0);
        fifun_ideal  = @(I0) spline(i_ideal, lambda, I0);

        % Coenergia
        cW_real   = integral(fifun_real,0,If);    
        cW_linear = integral(fifun_linear,0,If);
        cW_ideal  = integral(fifun_ideal,0,If);

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

    % cálculo da força
    F_reals   = diff(cW_reals)./diff(xs);
    F_linears = diff(cW_linears)./diff(xs);
    F_ideals  = diff(cW_ideals)./diff(xs);

    xfs = xs;
    xfs(1) = []; % Apagando o primeiro elemento para evitar problemas de convergência

    % Curva BxH
    figure('Name', 'Curva BxH do núcleo')
    plot(Hc_real,Bc,'x-',Hc_linear,Bc,'x-',Hc_ideal,Bc,'x-','LineWidth',1),zoom
    xlabel('H [ A/m ]','fontweight','bold','Fontsize',12)
    ylabel('B [ Wb/m² ]','fontweight','bold','Fontsize',12)
    legend('Curva Real','Curva Linear','Curva Ideal','Fontsize',16, location = 'southeast')
    axis([-50 500 0 2])
    grid 
    %set(gca, 'XScale', 'log')

%     % Curva λxi
%     ire = i_reals(:,1)';
%     figure('Name', 'Curva λxi')
%     plot(ire,lambda,'LineWidth',3),zoom
%     xlabel('i [ A ]','fontweight','bold','Fontsize',12)
%     ylabel("λ [ Wbe ]",'fontweight','bold','Fontsize',12)
%     grid minor

    % Deslocamento vs Coenergia
    figure('Name', 'Coenergia vs Deslocamento do êmbolo')
    plot(xs.*100,cW_reals,xs.*100,cW_linears,xs.*100,cW_ideals,'LineWidth',2),zoom
    xlabel('X [ cm ]','fontweight','bold','Fontsize',12)
    ylabel("W' [ J ]",'fontweight','bold','Fontsize',12)
    legend('Coenergia Real','Coenergia Linear','Coenergia Ideal','Fontsize',16)
    grid minor

    % Deslocamento vs Força
    figure('Name', 'Força vs Deslocamento do êmbolo')
    title('Deslocamento vs Força','Fontsize',16)
    plot(xfs.*100,F_reals,'-',xfs.*100,F_linears,'--',xfs.*100,F_ideals,'--','LineWidth',2),zoom
    xlabel('X [ cm ]','fontweight','bold','Fontsize',12)
    ylabel('F [ N ]','fontweight','bold','Fontsize',12)
    legend('Força Real','Força Linear','Força Ideal','Fontsize',16)
    grid minor

    prompt = input('Digite um novo valor de corrente ou 0 para encerrar o programa: '); % Calcula nova corrente ou encerra
    If = prompt;                % Corrente de excitação
    I  = linspace(0,If,n);      % Array de corrente para integração
    
    close all
end
