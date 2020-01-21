function [alphaArray, alphaArrayErrors] = calculate_alpha( RexArray, RexArray_E, dwXArray, B0Array, TempArray, QuantumCoherenceStringArray, OUTPUT_DEBUG_FLAG, varargin )
% Calculate exchange timescale alpha for each dispersion curve at each temperature
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%
% 2010/02/10 Create code
% 2011/04/19 Updated for GUARDD and classes
% 2011/06/30 Use dwX instead of B0
% 
% FUNCTION
%  Calculate exchange timescale alpha for each dispersion curve at each temperature
%  From Millet, O; Loria, PJ; Kroenke, CD; Pons, M; Palmer, AG III. 
%  "The Static Magnteic Field Dependence of Chemical Exchange Linebroadening
%   Defines the NMR Chemical Shift Time Scale." JACS, 122:2867-2877 (2000).
%
%  The exchange timescale can be helpful for determining the RD fitting equation required
% 
% INPUT VARIABLES
%  Rex
%  Rex_E
%  dwX array
%  B0
%  Temp
%  QuantumCoherenceString
%  OUTPUT_DEBUG_FLAG
%  varargin{1} = name of the title of debugging plot
%   
% OUTUPT VARS
%  Eq(13) and Eq(15)
%  0 <= alpha < 1  slow exchange
%  alpha = 1       intermediate exchange
%  1 < alpha <= 2  fast exchange
%  
% TO DO
%   Nothing!

if( ~iscell(QuantumCoherenceStringArray) )
    error('QuantumCoherenceStringArray must be a cell array of strings')
end

% Set up output plot if desired
if( OUTPUT_DEBUG_FLAG )    
    figure
    h = axes;
    hold(h, 'all');
end

% Initialize fitting variable arrays
Nc                      = length(RexArray);
alphaArray(1:Nc)        = NaN;
alphaArrayErrors(1:Nc)      = NaN;
interceptArray(1:Nc)    = NaN;

% Calculate alpha at each temperature / quantum coherence, if possible

% Itemize the unique temperatures (non-zero)
Temps   = unique(TempArray);
Temps   = Temps(Temps>0);
for t = 1:length(Temps)
    % Find all the curves at this temperature
    Ct = TempArray == Temps(t);
    
    if( OUTPUT_DEBUG_FLAG )    
        fprintf('\nWorking on Temp = %0.1f, curves: %s', Temps(t), sprintf('%d ',find(Ct)))    
    end
    
    % Itemize the unique quantum coherences at this temperature
    QuantumCoherences = unique(QuantumCoherenceStringArray(Ct));
    for q = 1:length(QuantumCoherences)        
        % Find all the curves at this quantum coherence
        Cq = strcmp(QuantumCoherenceStringArray, QuantumCoherences{q});
        
        % Find the curves that are at this Temperature AND QuantumCoherence
        Ctq = and(Ct, Cq);
        
        if( OUTPUT_DEBUG_FLAG )    
            fprintf('\n\tWorking on Quantum Coherence = %s, curves: %s', ...
                QuantumCoherences{q}, sprintf('%d ',find(Cq)))
            fprintf('\n\t\tCurves at uniform Temp and QC: %s', sprintf('%d ',find(Ctq)))
        end        
        
        %% Extract the data for fitting
        % Calculate natural log to fit data to line
        %  Error propagation: E{ a*ln(k) } = a*E{k}/k
        %X   = log(B0Array(Ctq));
        X   = log(dwXArray(Ctq));
        Y   = log(RexArray(Ctq));
        Y_E = RexArray_E(Ctq) ./ RexArray(Ctq);
        
        % Remove +/-Inf data points (Rex=0 -> log(0)=-Inf)
        iInf = isinf(Y);
        if( any(iInf) )                            
            X   = X(~iInf);
            Y   = Y(~iInf);
            Y_E = Y_E(~iInf);
        end
        
        % Fit the data
        [ m, m_E, b, b_E ] = fit_line( X, Y, Y_E );

        % Save the result for each curve at this temperature / QC
        alphaArray(Ctq)        = m;
        alphaArrayErrors(Ctq)  = m_E;
        interceptArray(Ctq)    = b;
        
        % Output fit result if desired
        if( OUTPUT_DEBUG_FLAG )
            iCtq_nonzero = find(Ctq);
            index = iCtq_nonzero(1);
            
            errorbar( X, Y, Y_E, 'o', 'MarkerSize', 10, 'Linewidth', 2 );
            text( X(1), Y(1), sprintf('     \\alpha=%0.3f\\pm%0.3f at %d^oC (%s)', ...
                alphaArray(index), alphaArrayErrors(index), ...
                Temps(t)-273, QuantumCoherences{q}), 'FontSize', 16);
        end
    end
end

%% Plot the fit on top of data
if( OUTPUT_DEBUG_FLAG )
    XLIM = get(h, 'XLim');
    YLIM = get(h, 'YLim');
    
    for t = 1:length(Temps)
        % Find all the curves at this temperature
        Ct = TempArray == Temps(t);

        % Itemize the unique quantum coherences at this temperature
        QuantumCoherences = unique(QuantumCoherenceStringArray(Ct));
        for q = 1:length(QuantumCoherences)        
            % Find all the curves at this temperature and quantum coherence
            %Ctq = find(strcmp(QuantumCoherenceStringArray(Ct), QuantumCoherences(q)));
            Cq = strcmp(QuantumCoherenceStringArray, QuantumCoherences{q});
            
            % Find the curves that are at this Temperature AND QuantumCoherence
            Ctq = and(Ct, Cq);
            
            % Get the first curve for plotting
            iCtq_nonzero = find(Ctq);
            index = iCtq_nonzero(1);

            % Plot the linear fit using the first curve C(1)
            X = [-5, 10];
            plot( X, X.*alphaArray(index) + interceptArray(index), '--k');
        end
    end

    set(h,'XLim',XLIM);
    set(h,'YLim',YLIM);

    %grid(h, 'on');
    box(h,'on');

    set(h, 'FontSize', 16)
    %xlabel(h, 'Ln( B_0 (MHz) )', 'FontSize', 14)
    xlabel(h, 'Ln( \Delta\omega_X (Hz) )', 'FontSize', 16)
    ylabel(h, 'Ln( R_{Ex} (Hz) )', 'FontSize', 16)
    
    if( nargin == 8 )
        title(h,varargin{1}, 'FontSize', 16);
    end
end
