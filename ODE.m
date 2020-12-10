%% Deterministic Modeling:
% Plotting ODE of mRNA Counts:
syms x
%Rates (per minute, Given by Professor Golding)
kR = 10; % transcription rate mRNA
degR = 1/3;% degradation rate mRNA
kP = 6; % transcription rate protein
degP = 1/60; % degradation rate protein
a = (kP*kR)/degR; % expression rate
% Equation
figure()
fplot(-1*kR*(1/degR)*exp(-degR*x) + (1/degR)*kR)
hold on
fplot(-1*(a/degP)*exp(-degP*x) + a/degP)
% Steady State
fplot(kR/degR,'--','LineWidth', 2)
fplot(a/degP,'--','LineWidth', 2)

axis([0 1000 0 12000])
legend({'no. of mRNA','Steady State Concentration'},'Location','northwest')
xlabel('time(minutes)')
ylabel('no. of mRNA')
title('Deterministic modelling of mRNA counts')


