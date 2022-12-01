SMSE = readmatrix('SMSE_estimate.txt');
RMSE = readmatrix('RMSE_estimate.txt');
Rho = readmatrix('Rho_estimate.txt');
pRho = readmatrix('pRho_estimate.txt');
MSLL = readmatrix('MSLL_estimate.txt');
EXPV = readmatrix('EXPV_estimate.txt');

NMacc = table(SMSE, RMSE, Rho, pRho, MSLL, EXPV);