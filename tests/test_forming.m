% Copyright 2019, Eduardo L. T. Conceicao
% Available under the GPL-3


rng('default');

%% Test massive material
randnet = forming(40, double(intmax('int64')), 200, 200, nfib = 10000);
nopores = porosity(randnet);
assert(nopores.interfiber == 0)

%% Test hydrodynamic smoothing
randnet = forming(40, 1, 200, 200, ...
                  nfib = 10000, acceptanceprob = 0.1);
sheet_thicker = thickness(randnet);
randnet = forming(40, 1, 200, 200, ...
                  nfib = 10000, acceptanceprob = 0.01);
sheet_less_thick = thickness(randnet);
assert(sheet_thicker.apparent >= sheet_less_thick.apparent)

%% Test fiber flexibility
randnet = forming(40, 1, 200, 200, nfib = 10000);
sheet_thicker = thickness(randnet);
randnet = forming(40, 4, 200, 200, nfib = 10000);
sheet_less_thick = thickness(randnet);
assert(sheet_thicker.apparent >= sheet_less_thick.apparent)
