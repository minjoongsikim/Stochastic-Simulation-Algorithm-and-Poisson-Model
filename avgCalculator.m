function [ t, x, x2, avg1, avg2 ] = avgCalculator( stoich_matrix, propensity_fcn, tspan, x0,...
                                  params, output_fcn, MAX_OUTPUT_LENGTH)

if ~exist('MAX_OUTPUT_LENGTH','var')
    MAX_OUTPUT_LENGTH = 1000000;
end
if ~exist('output_fcn', 'var')
    output_fcn = [];
end
if ~exist('params', 'var')
    params = [];
end
    
%% Initialize
%num_rxns = size(stoich_matrix, 1);
num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;
count = 0;
countM = 0;
new_rxn = 0;
%% MAIN LOOP
while T(rxn_count) < tspan(2)        
    % Calculate reaction propensities
    a = propensity_fcn(X(rxn_count,:), params);
    
    % Sample earliest time-to-fire (tau)
    a0 = sum(a);
    r = rand(1,2);
    tau = -log(r(1))/a0; %(1/a0)*log(1/r(1));
    
    % Sample identity of earliest reaction channel to fire (mu)
    [~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]); 
    
    % ...alternatively...
    %mu = find((cumsum(a) >= r(2)*a0), 1,'first');
    
    % ...or...
    %mu=1; s=a(1); r0=r(2)*a0;
    %while s < r0
    %   mu = mu + 1;
    %   s = s + a(mu);
    %end
%     if rxn_count + 1 > MAX_OUTPUT_LENGTH
%         t = T(1:rxn_count);
%         x = X(1:rxn_count,:);
%         warning('SSA:ExceededCapacity',...
%                 'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
%         return;
%     end
    
    % Update time and carry out reaction mu
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
    rxn_count = rxn_count + 1;
    if T(rxn_count) > 400
        disp(T(rxn_count))
        count = count + X(rxn_count);
        countM = countM + X(rxn_count,2);
        new_rxn = new_rxn + 1;
    else 
        disp(T(rxn_count))        
    end
    
    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count,:)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
            warning('SSA:TerminalEvent',...
                    'Simulation was terminated by OutputFcn.');
            return;
        end 
    end
end  

   
% Return simulation time course
t = T(1:rxn_count);
x = X(1:rxn_count,1);
x2 = X(1:rxn_count, 2);
avg1 = count/new_rxn;
avg2 = countM/new_rxn;
fprintf('avg = %f\n', count) 
fprintf('total avg Protein = %f\n', count/new_rxn) 
fprintf('total avg mRNA = %f\n', countM/new_rxn)
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,1) = X(rxn_count-1,1);
    x2(end, 1) = X(rxn_count-1, 2);
    
end
end


