function [Kappa, Acc] = GetKappaAcc(tp2)

a = sum(tp2);
b = sum(tp2.');
N = sum(sum(tp2));
accSum = sum(diag(tp2));

p0 = accSum/N;
pe = a*b'/(N^2);

Kappa = (p0-pe)/(1-pe);
Acc = sum(diag(tp2))/sum(sum(tp2));