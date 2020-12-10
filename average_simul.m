function average_simul()
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
total1 = 0;
total2 = 0;
s = 10; %number of simuls for avg.
avgArray1 = zeros(s());
avgArray2 = zeros(s());
%% Specify reaction network
pfun = @propensities_2state;
stoich_matrix = [ 1  0    %transcription
                  0  1    %translation
                  -1  0    %mRNA decay
                  0  -1 ]; %protein decay

%% Run AVGsimulation
for i = 1:s
    [t,x,x2, avg1, avg2] = avgCalculator(stoich_matrix, pfun, tspan, x0, p);
    avgArray1(i) = avg1;
    avgArray2(i) = avg2;
    total1 = total1 + avg1;
    total2 = total2 + avg2;
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
fprintf('total 1 = %f\n', total1) 
fprintf('total 1 = %f\n', total2)
fprintf('total 1 = %f\n', total1) 
fprintf('array1 = %f\n', avgArray1)
fprintf('array2 = %f\n', avgArray2) 
std1 = std(avgArray1);
std2 = std(avgArray2);
fprintf('std2 = %f\n', std1)
fprintf('std2 = %f\n', std2) 

end