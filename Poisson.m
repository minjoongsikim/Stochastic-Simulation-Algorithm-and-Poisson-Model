
% Use Distributed Data.
k = poissrnd(60,1000,1);
% Fit a poisson distribution
pd = fitdist(k,'poisson');
% Plot histogram and fit
figure
hold on
histogram(k,0:100,'Normalization','probability')

x = 0:100;
lambda = 60;
f = (lambda.^x) .* exp(-lambda) ./ factorial(x);
plot(x,f, "-", 'LineWidth', 2)
legend('Distribution of mRNA', 'Poisson Distribution \lambda = 60')
title('Poisson Distribution of steady-state mRNA values')
xlabel('number of mRNA') 
ylabel('Probability') 