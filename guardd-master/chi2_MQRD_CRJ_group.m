function [chi2_tot, df, chi2Array] = chi2_MQRD_CRJ_group( fitResult, p_fmincon, varargin )
% Calculate chi^2 for fit to entire group
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2008/07/07 Start code
% 2008/08/14 Preallocate arrays for speed
% 2008/09/05 Changed vcpmg_matrix so that values are row vectors, not columns
% 2008/09/24 Now returns array of SSE for each fit
% 2009/02/09 Returns chi2 = SUM( (y-fit)^2 / (ey)^2 )
%            This is so that fit discrepancies are normalized to data's error
% 2009/03/09 Accepts p_id_matrix and Cp_id_matrix instead of C_matrix
% 2011/01/12 Vectorized some code (for speed)
% 2011/04/18 Use group class structure
% 2011/04/30 Constrain temperature-depdendence of rates CONSTRAIN_RATE_ANALYSIS
% 2011/06/06 Use fitResult instead of Group for optional CONSTRAIN_RATE_ANALYSIS
% 2012/01/11 Return chi2Array -- the chi2 value for each curve
%
% This will return the sse from a global fit of the MRD-MQ data
% This is the multi-dimensional analog of mse_fit_mq_r2eff
%
%
% INPUT VARIABLES
%   fitResult -> Contains info on group(data and constraints)
%   p     -> linearized parameter array for fitting (must be unpacked to
%            matrix form via group)
%   R2eff_Surrogate{ctot} -> Contains alternate Y data for each curve in group
%                       Used for Monte Carlo error estimation 
%                       Set as [] if not needed (or do not supply argument)
%
% RETURN VARIABLES
%  chi2_tot     = Global chi^2 for fitting all curves
%  df           -> Degrees of freedom in the fit
%
OUTPUT_DEBUG = false;    

% Processing optional argument for R2eff_Surrogate
if( nargin == 3 )
    R2eff_Surrogate = varargin{1};                
else
    R2eff_Surrogate = [];
end

group   = fitResult.parentGroup;
Nctot   = group.getNumCurves();
Np      = 5;

% De-linearize parameter array to matrix form
pMatrix = group.delinearizePFmincon(p_fmincon, fitResult);

if( OUTPUT_DEBUG )
    clc
    p_fmincon(1)
    p_fmincon(2)
    pMatrix
end

% Iterate through each curve in the group
chi2Array = zeros(1,Nctot);
NobsTotal = 0;
for ctot = 1:Nctot
    % Get curveset and curve number
    [cs,c] = group.getCurvesetCurve( ctot );
    
    % Fit each curve
    vcpmg   = group.curvesets{cs}.curves{c}.vcpmg;
    
    % Get the alternate Y-axis data if available
    if( ~isempty(R2eff_Surrogate) )
        R2eff = R2eff_Surrogate{ctot};
        if( length(R2eff) ~= length(vcpmg) )
            error('Supplied R2eff (Y) data does not match the size of vcpmg (X) data');
        end        
        
    % Otherwise, get the data from the group itself
    else
        R2eff   = group.curvesets{cs}.curves{c}.R2eff;
    end
    
    eR2eff  = group.curvesets{cs}.curves{c}.eR2eff;
    TCPMG   = group.curvesets{cs}.curves{c}.TCPMG;
    
    params  = pMatrix(ctot,1:Np);
    Nfixed  = 0;
    
    %{
    vcpmg'
    R2eff'
    eR2eff'
    TCPMG
    params
    Nfixed
    
    chi2_MQRD_CRJ( vcpmg, R2eff, eR2eff, TCPMG, params, Nfixed )
    input('\nHit return');
    %}
    
    [chi2Array(ctot), void, void]   = chi2_MQRD_CRJ( vcpmg, R2eff, eR2eff, TCPMG, params, Nfixed );
    NobsTotal                       = NobsTotal + group.curvesets{cs}.curves{c}.Nobs;
    
    if( OUTPUT_DEBUG )
        fprintf('\n\nFitted curve %d\tchi2 = %0.1f', ctot, chi2Array(ctot));
        fprintf('\n\tp\t= %s',sprintf('%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f', params))
        

        %fprintf('\n\teR2eff\t= %s\n\tCHI2\t= %s', ...
        %    sprintf('%2.2f\t',eR2eff_matrix(c,:)), sprintf('%2.2f\t',chi2(c)));
        input '\nEnter to continue chi2_MQRD_CRJ_group';
    end
end

% Total chi2 for all the curves
chi2_tot = sum(chi2Array);

% Degrees of freedom in global fit = Nobservations - Nparams
df = NobsTotal - length(p_fmincon);

if( OUTPUT_DEBUG )    
    fprintf('\n\nDF=%2.2f\tCHI2_total=%f', df, chi2_tot);
    input '\nDone with chi2_MQRD_CRJ_group'
end
