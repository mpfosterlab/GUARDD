function [ m, m_E, b, b_E ] = fit_line( X, Y, Y_E )
%fit_line Fits a line Y = m*X + b to X,Y data with errors Y_E
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2011/05/?? Started code
% 2011/05/25 Updated to get unique observations only

% Make sure data are column vectors
if( size(X,2) > 1 )
    X   = X';
    Y   = Y';
    Y_E = Y_E';
end

% Get unique observations only
[VOID, i_unique, VOID] = unique(X);
X   = X(i_unique);
Y   = Y(i_unique);
Y_E = Y_E(i_unique);

% Number of observations
Nobs = length( unique(X) );

% Must have more than one unique point to fit the line
if( Nobs > 1 )
    % If there are no errors, or any errors are NaN set weights all to unity
    if( any(Y_E==0) || any(isnan(Y_E)) || any(isinf(Y_E)) || Nobs==2 )
        weights = ones(Nobs,1);
    else                
        weights = 1 ./ Y_E;
    end

    s = fitoptions( 'Method','NonlinearLeastSquares',...
                    'Weights', weights, ...
                    'Lower',[-Inf,-Inf], ...
                    'Upper',[Inf,Inf], ...
                    'Startpoint',[1 1] );
    fit_type = fittype('m*x + b','options',s);            
    [params, outputfit] = fit(X, Y, fit_type);

    %% Estimate error in fit
    % Only get confidence interval if there are more than 2 points for line
    if( Nobs > 2 )
        outputfit_conf      = confint(params, 0.683); % Get limits at 1 sigma
        b1                  = outputfit_conf(1);    % Lower limit at 1 sigma
        b2                  = outputfit_conf(2);    % Upper limit at 1 sigma
        b                   = 0.5*(b1+b2);          % Fitted value
        b_E                 = abs(b2-b);            % Error in fit (1 sigma)
        m1                  = outputfit_conf(3);    % Lower limit at 1 sigma
        m2                  = outputfit_conf(4);    % Upper limit at 1 sigma
        m                   = 0.5*(m1+m2);          % Fitted value
        m_E                 = abs(m2-m);            % Error in fit (1 sigma)

    % Otherwise construct a synthetic error from error in observed variable
    else
        synth_f_error       = sqrt(sum( (Y_E ./ Y).^2) );
        b                   = params.b;
        b_E                 = b*synth_f_error;
        m                   = params.m;
        m_E                 = m*synth_f_error;
    end
    
% Not enough points to fit the line
else
    m   = NaN;
    m_E = NaN;
    b   = NaN;
    b_E = NaN;
end

