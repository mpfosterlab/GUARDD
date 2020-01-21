function Rex = calculate_Rex( Tcpmg, p_sim )
% Calculate exchange broadening Rex for the dispersion
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/01/28 Create code
% 2011/04/22 Ensure Rex >= 0 (round 1e-3 to zero)
% 2017/05/26 Fixed a bug with not properly dealing with R2=Inf
% 
% FUNCTION
%  Calculate exchange broadening Rex for the dispersion
%  Rex = R2Eff(vCPMG=Inf) - R2Eff(vCPMG=0) ~= R2Eff(vCPMG = 10kHz) - R2Eff(vCPMG = 0.1 Hz)
% 
% INPUT VARS
%  Tcpmg      = Total CPMG time for the dispersion curve (sec)
%  p_sim      = Parameters for simulating single dispersion curve
%   
% OUTUPT VARS
%  Rex
%  
% TO DO
%   Nothing!

Np  = 5;

% Try to calculate at vcpmg = 0.1. If this does not work (NaN) add 0.1 and try again
R20     = NaN;
vcpmg0  = 0;
while( isnan(R20) || isinf(R20) )
    vcpmg0  = vcpmg0 + 0.1;
    R20     = model_MQRD_CRJ(vcpmg0, Tcpmg, p_sim(1:Np));
end

% Try calculate at vcpmg = 10000, if this does not work, reduce by 1 and
% try again
R2inf       = NaN;
vcpmginf    = 10000;
while( isnan(R2inf) || isinf(R2inf) )
    vcpmginf    = vcpmginf - 1;
    R2inf       = model_MQRD_CRJ(vcpmginf, Tcpmg, p_sim(1:Np));
end

% Done!
Rex = R20 - R2inf;

if( Rex < 1e-3 )
    Rex = 0;
end