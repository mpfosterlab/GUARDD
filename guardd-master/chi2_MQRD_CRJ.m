function [chi2, chi2red, e] = chi2_MQRD_CRJ( vcpmg, R2eff, eR2eff, TCPMG, p, Nfixed )
% Calculate reduced chi^2 between observed and simulated MQ-CPMG RD data
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2008/05/30 Start code
% 2008/08/15  Vectorized calculations of all arrays for speed
% 2008/09/05 Changed vcpmg so it is a row vector (not a column vector)
% 2011/04/14 Rename to chi2_MQ_CRJ_RD
% 2011/06/07 Prevent divide by zero if eR2eff = 0, instead assume eR2eff=1
%            (i.e., chi2 is actually sum of squared errors)
%
% PURPOSE
%  Calculate reduced chi^2 between observed and simulated MQ-CPMG RD data
%
%  Equation numbers refer to equations from
%   D.M. Korzhnev, K. Kloiber, V. Kanelis, V. Tugarinov, L.E. Kay, Probing slow dynamics in high
%   molecular weight proteins by methyl-trosy nmr spectroscopy: application to a 723-residue enzyme,
%   J. Am. Chem. Soc. (2004) 3964-3973.
%
% INPUT VARS
%  Measured and/or known data
%   vcpmg   = CPMG frequency array
%   R2eff   = R2Eff relaxation rate array
%   TCPMG   = Total CPMG time
%   p       = Parameter array (dwH,dwX,Pa,Kex,R2mq)
%   Nfixed    Number of parameters not varied in the fit (to calculate df -> reduced chi2)
%
% OUTPUT VARS
%  Note: mse and sse are calculated to be reduced chi2 and chi2, respectively
%  If the function is called with 1 return value, it will be chi2
%  If the function is called with 3 return values, it will be [chi2, chi2red, e]
%
% LIMITATIONS
%  This eqn can't be used for SQ RD unless TCPMG > ~30ms (Korzhnev, et al., 2004)

% Number of observations
N = length(vcpmg);
%vcpmg = vcpmg'; 09/05/08

% Set values of the parameters to vary for best fit
% If five parameters are to be fit
dw_h	= p(1);
dw_x	= p(2);
pa      = p(3);
pb      = 1 - p(3);
kex     = p(4);
R2mq	= p(5);

% Degrees of freedom
% = (Num of obs) - (Num of model parameters) + (Num parameters fixed)
df = N - length(p) + Nfixed;

% These arrays accompany R2eff values
delta = 1 ./ (4*vcpmg);
n = TCPMG ./ (4*delta);

% (3.10) scalars
z_plus  = dw_h - dw_x + sqrt(-1)*kex;
z_minus = dw_h - dw_x - sqrt(-1)*kex;

d_plus  = dw_h + dw_x + sqrt(-1)*kex;
d_minus = dw_h + dw_x - sqrt(-1)*kex;


% (3.9) vector
mZ = -sqrt(-1)*kex*sqrt(pa*pb)*( d_minus - 2*dw_x*sin(d_minus.*delta) ./ ...
    sin((d_minus+z_minus).*delta) )./(d_minus*z_minus);

% (3.8) vector
mD = sqrt(-1)*kex*sqrt(pa*pb)*( z_plus + 2*dw_x*sin(z_plus.*delta) ./ ...
    sin((d_plus+z_plus).*delta) ) ./ (d_plus*z_plus);

% (3.7) - vector
Q = real(1 - mD.^2 + mD.*mZ - mZ.^2 + 0.5*(mD + mZ).*sqrt(pb/pa));

% (3.6) - scalar
zeta = -2*dw_x*(sqrt(-1)*dw_h + (pa-pb)*kex);

% (3.5) - scalar
Psi = (sqrt(-1)*dw_h + (pa-pb)*kex)^2 - dw_x^2 + 4*pa*pb*kex*kex;

% (3.4) - vectors
eta_plus = sqrt(2)*delta*sqrt( sqrt(Psi^2 + zeta^2) + Psi );
eta_minus= sqrt(2)*delta*sqrt( sqrt(Psi^2 + zeta^2) - Psi );

% (3.3) - scalars
D_minus= 0.5*( (Psi + 2*dw_x^2) / sqrt(Psi^2 + zeta^2) - 1 );
D_plus = 0.5*( (Psi + 2*dw_x^2) / sqrt(Psi^2 + zeta^2) + 1 );


% (3.2) - vector
lambda1 = R2mq + 0.5*(kex - (1./(2*delta)) .* ...
    acosh( D_plus*cosh(eta_plus) - D_minus*cos(eta_minus) ) );

% (3.1) - vector
mR2eff = real(lambda1) - log(Q)./(4*n.*delta);

% Get the residuals
e   = R2eff - mR2eff;

sse = sum(e.^2);
mse = sse / df;

% Prevent divide by zero while computing chi2
% Instead, assume eR2eff=1 (i.e., chi2 is actually sum of squared errors)
if( any(eR2eff == 0) )
    chi2    = sum( (e ./ 1).^2 );
else
    chi2    = sum( (e ./ eR2eff).^2 );
end
chi2red = chi2 / df;

%    fprintf(1, '\nsse=%f\tdw_h=%f\tp(1)=%f', sse, dw_h, p(1));    
%    fprintf(1, '\ne = %f\tsse = %f\tmse = %f',e,sse, mse);
