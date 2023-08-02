% Tabela e funções - Questão 3

% Pontos da curva BxH real
Hreal = 1000.*[-50.5, -50, -48.5, -47.5, -45.5, -40, -30, -20, -10, 0];
Breal = [0, 0.1, 0.4, 0.6, 0.8, 0.98, 1.08, 1.13, 1.17, 1.2];

% Pontos da curva BxH linear
Hlinear = [Hreal(1),Hreal(length(Hreal))];
Blinear = [Breal(1),Breal(length(Breal))];

% Funções BxH
Hfun_real = @(B) spline(Breal, Hreal, B);
Hfun_linear = @(B) spline(Blinear, Hlinear, B);

% Ponto de coercitividade
Hcoerc = Hreal(1);

% Pontos de operação normal
Bopn = 1.07;
Hopn_real = Hfun_real(Bopn);
Hopn_linear = Hfun_linear(Bopn);

clear Hreal Breal Hlinear Blinear

% n = 200;
% Bc = linspace(0,1.2,n);
% figure('Name', 'Curva BxH do ímã')
% title('Curva BxH do ímã')
% plot(Hfun_real(Bc),Bc,'b-',Hfun_linear(Bc),Bc,'r--')
% axis([-55000 0, 0 1.2])
% grid minor