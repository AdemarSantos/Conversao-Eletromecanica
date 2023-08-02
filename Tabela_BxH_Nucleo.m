% Tabela e funções - Questões 1 e 2

% Pontos da curva BxH real
Hreal = [0, 68, 135, 203, 271, 338,...
        406, 474, 542, 609, 1100,...
        1500, 2500, 4000, 5000, 9000,...
        12000, 20000, 25000];

Breal = [0.000, 0.733, 1.205, 1.424, 1.517,...
        1.560, 1.588, 1.617, 1.631, 1.647,...
        1.689, 1.703, 1.724, 1.731, 1.738,...
        1.761, 1.770, 1.800, 1.816];

% Pontos da curva BxH linear    
Hlinear = [0, 68];
Blinear = [0.000, 0.733];

% Funções BxH
Hfun_real = @(B) spline(Breal, Hreal, B);
Hfun_linear = @(B) spline(Blinear, Hlinear, B);

clear Hreal Breal Hlinear Blinear
