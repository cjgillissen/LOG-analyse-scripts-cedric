%% Noise correlations for widefield data

L1 = length(figcorr);
L2 = length(figerr);
labels = [ones(L1,1);zeros(L2,1)];
scores = [figcorr;figerr];
[X,Y,T,AUC1] = perfcurve(labels,scores,1);
 
%Assuming figcorr and figerr are vectors of the responses on correct and error trials
%respectively.
 
Cheers, Matt