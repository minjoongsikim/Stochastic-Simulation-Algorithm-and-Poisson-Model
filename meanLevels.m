
title('Avg Values of mRNA and Protein')
hold on
yyaxis left
ylim([0 40])
X = categorical({'mRNA Steady State', 'mRNA Gillespie Avg.'});
%% ADD ERROR BARS
vals = [30 30.278];
bar(X,vals)
err1 = [.607 0];
err3 = [.607 0];
er = errorbar(1,vals(2),err1(1));    
er.Color = [0 0 0];  
er.CapSize = 15;
er.LineStyle = 'none';



yyaxis right
ylim([0 12000])
Z = categorical({'protein Steady State', 'protein Gillespie Avg.'});
vales = [10800 10773.387];
bar(Z,vales)

%% ADD ERROR BARS
q = 3:4;
err2 = [184.25 0];
err0 = [184.25 0];
er2 = errorbar(3, vales(2),err2(1));
er2.CapSize = 15;
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  
