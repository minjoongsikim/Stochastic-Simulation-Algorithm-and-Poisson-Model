function ssa_100_simul()
% Simulate a two-state model of gene expression
import Gillespie.*
warning('on')
%% Reaction network:
%   1. transcription:       0       --kR--> mRNA
%   2. translation:         mRNA    --kP--> mRNA + protein
%   3. mRNA decay:          mRNA    --gR--> 0
%   4. protein decay:       protein --gP--> 0
%% Rate constants
p.kR = 10; % rate of transcription mRNA      
p.kP = 6;  % rate of translation protein                    
p.gR = 1/3; % rate of degradation mRNA
p.gP = 1/60; % rate of degradation protein
%% Rates (per minute, Given by Professor Golding)
kR = 10; % transcription rate mRNA
degR = 1/3;% degradation rate mRNA
kP = 6; % transcription rate protein
degP = 1/60; % degradation rate protein
a = (kP*kR)/degR; % expression rate
%% Initial state
tspan = [0, 1000]; %seconds
x0    = [0, 0];     %mRNA, protein
%% Specify reaction network
pfun = @propensities_2state;
stoich_matrix = [ 1  0    %transcription
                  0  1    %translation
                  -1  0    %mRNA decay
                  0  -1 ]; %protein decay

%% Run simulation
[t,x,x2] = directMethod100(stoich_matrix, pfun, tspan, x0, p);
%% Plot time course
figure();
hold on
time = 0:1000; % time range
% %% Plot Deterministic Equations
% w = -1*kR*(1/degR)*exp(-degR*time) + (1/degR)*kR; % mRNA deterministic equation
% yyaxis left
% ylim([0, 100])
% plot(time,w,'k--','linewidth',1)
% r = -1*(a/degP)*exp(-degP*time) + a/degP;
% yyaxis right
% ylim([0, 15000])
% plot(time,r,'--')
%% Plot Steady States
steadyR = kR/degR; % Steady State mRNA
yyaxis left
plot(time,steadyR,'r:')
yyaxis right
steadyP = a/degP; % Steady State Protein
plot(time, steadyP,'k:')
%% Plot Stochastic Models
yyaxis left
ylabel('number of mRNA');
h1 = stairs(t,x, 'b:','linewidth',.1); set(gca,'XLim',tspan);
yyaxis right
h2 = stairs(t,x2, ':','linewidth', 1); set(gca,'XLim',tspan);
xlabel('time (s)');
ylabel('number of protein');
%% Plot Deterministic Models

w = -1*kR*(1/degR)*exp(-degR*time) + (1/degR)*kR; % mRNA deterministic equation
yyaxis left
ylim([0, 100])
plot(time,w,'k--','linewidth',.01);
r = -1*(a/degP)*exp(-degP*time) + a/degP; %Protein deterministic equation
yyaxis right
ylim([0, 15000])
plot(time,r,'k--', 'linewidth', 1);
lgd = legend([h1 h2],'mRNA', 'protein');
lgd.FontSize = 14;
title('100 Stochastic Simulations and Differential Model')
end
function a = propensities_2state(x, p)
% Return reaction propensities given current state x
mRNA    = x(1);
protein = x(2);
a = [p.kR;            %transcription
     p.kP*mRNA;       %translation
     p.gR*mRNA;       %mRNA decay
     p.gP*protein];   %protein decay
end