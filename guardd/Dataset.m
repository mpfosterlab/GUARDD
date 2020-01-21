classdef Dataset < handle
    % Dataset class: for a single RD dataset
    %
    % (C) Ian Kleckner [ian.kleckner@gmail.com]
    %  Foster Lab, The Ohio State University
    % GUARDD software [http://code.google.com/p/guardd/]
    %  GNU GPL3 License
    %
    % 2011/01/17 Start coding
    % 2011/06/16 Fixed error calculating eIntensity from eR2eff in addData()
    % 2017/04/25 Added 19F (IK)
    %
    % Properties and methods for a single RD datset
    
    properties (SetAccess = private)
        % SetAccess = private properties can be read via ".", but not set
        
        name        = 'DatasetName';      % Name of dataset (e.g., 'MQ 800MHz 25C')
        AX_String   = 'AX';     % 13C or 15N
        B0          = 0;        % Magnetic field strength (MHz)
        Temp        = 0;        % Temperature (K)
        Temp_C      = 0;        % Temperature (C)
        TCPMG       = 0;        % Total CPMG time in pulse sequence (sec)
        SQX         = false;    % (boolean) data set is SQ for nucleus X

        Intensity_matrix = [];       % 2D array of signal intensity values (-)
        eIntensity_matrix= [];       % 2D array of errors in signal intensity (-)
        R2eff_matrix     = [];       % 2D array of R2eff values (Hz)
        eR2eff_matrix    = [];       % 2D array of errors in R2eff (Hz)
        vcpmg            = [];       % Array of vcpmg values (Hz)        
        Nobs             = 0;        % Number of observations (length(vcpmg))
        Nc               = 0;        % Number of curves in the data set
        
        index_array     = {};   % Array of residue numbers
        atom_array      = {};   % Array of names of atoms (N, H\alpha, C\delta_1, etc.)
        residue_array   = {};   % Array of names of residues (Ile, Leu, Arg, etc.)
        
        nlin_filename   = '';   % Name of nlin.tab file for input data (NMRPipe)
        vcpmg_filename  = '';   % Name of file for vcpmg values
        %r2table_filename= '';   % Name of r2table formatted file (R2eff, and vcpmg)
        
        curves          = {};
        
        % Cell matrix for parameters to describe this class via table
        %  Used in generateParamTable()
        param_info = {...
            'name:(String) Name of dataset (e.g., MQ 800MHz 25C)', ...
            'AX_String:(String) 13C, 15N, or 19F', ...
            'B0:(MHZ) Magnetic field strength', ...
            'Temp_C:(C) Temperature', ...
            'TCPMG:(sec) Total CPMG time in pulse sequence', ...
            'SQX:(Boolean) Data set is SQ for nucleus X' ...
                      };
                  
                  
        parentSession   = [];   % The session that owns this dataset
    end
    
    methods
        function obj = Dataset( parentSession )
            %% Default constructor
            obj.parentSession = parentSession;
        end
        
        function addCurve( obj, curve )
           %% Add a curve to the dataset (e.g., from an individual file) 
            
        end
        
        function addData(obj, INDEX, ATOM_STRING, RESIDUE_STRING, X, Y, Y_E, Y_DATA_STRING )
            %% Add a single RD curve to the dataset
            % X,Y,Y_E are vectors for data (multiple observations)
            % Y_E = [] => Use multiple observations for calculating errors
            % Y_E = (x x x x) => Specify errors explicitly
            % Note: Use Y_E = fError .* Y for fractional error (fError is a number)
            % 2011/06/09 Add RESIDUE_STRING

            % Enable output for debugging purposes
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_ADD_DATA;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            % Check for valid data type
            switch upper(Y_DATA_STRING)
               case {'R2', 'INTENSITY'}

               otherwise
                   error('Invalid Y_DATA_STRING specified, "%s". See code for options', Y_DATA_STRING);
            end   

            % Calculate errorrs from duplicate observations
            if( isempty(Y_E) )
               [X, Y, Y_E] = Dataset.calculateErrorsUsingDuplicates( X, Y );
               
               if( OUTPUT_DEBUG )
                   fprintf('\nCalcualated errors');
                   fprintf('\n\tX\tY\tY_E\n');
                   disp( [X' Y' Y_E'] )
               end
            end

            % Assign the data based on its type
            switch upper(Y_DATA_STRING)
               case 'R2'
                   VCPMG    = X;
                   R2eff    = Y;
                   eR2eff   = Y_E;

                   % Calculate Intensity from R2Eff
                   if( obj.TCPMG > 0 )
                       I0           = 1;
                       Intensity    = I0 * exp( -obj.TCPMG .* R2eff );                       
                       eIntensity   = I0 * exp(-obj.TCPMG .* R2eff) * obj.TCPMG .* eR2eff;
                   else
                       error('R2eff requires defined Intensity and TCPMG'); 
                   end

               case 'INTENSITY'
                   VCPMG        = X;
                   Intensity    = Y;
                   eIntensity   = Y_E;
                   
                   % Find/remove VCPMG = 0 points
                   oZero = find(VCPMG == 0);
                   if( ~isempty(oZero) )
                       % Normalize each observation
                       I0           = Intensity(oZero);
                       Intensity    = Intensity ./ I0;
                       eIntensity   = eIntensity ./ I0;
                       
                       if( OUTPUT_DEBUG )
                           fprintf('\nFound I0 = %0.1f at observation %d', I0, oZero);
                       end                       
                   end
                   
                   % Only keep the non-zero VCPMG observations
                   oNonZeros    = find(VCPMG ~= 0);
                   VCPMG        = VCPMG(oNonZeros);
                   Intensity    = Intensity(oNonZeros);
                   eIntensity   = eIntensity(oNonZeros);

                   % Calculate R2eff from intensity
                   if( obj.TCPMG > 0 )
                       R2eff        = -log( Intensity ) ./ obj.TCPMG;
                       eR2eff       = eIntensity ./ ( Intensity * obj.TCPMG );
                   else
                       error('R2eff requires defined Intensity and TCPMG'); 
                   end
            end

            % Count the new curve
            Nc = obj.Nc + 1;

            % Number of observations
            Nobs = length(VCPMG);

            % If this is the first curve in the dataset, then set the
            % observation array
            if( Nc == 1 )
               obj.Nobs     = Nobs;           
               obj.vcpmg    = VCPMG;

            elseif( Nobs ~= obj.Nobs )
               error('All curves in a single dataset must have same number of observations (VCPMG values)');               
            end
            
            % Add these data to the data matrices
            obj.Intensity_matrix(Nc,1:Nobs)     = Intensity;
            obj.eIntensity_matrix(Nc,1:Nobs)    = eIntensity;
            obj.R2eff_matrix(Nc,1:Nobs)         = R2eff;
            obj.eR2eff_matrix(Nc,1:Nobs)        = eR2eff;

            % Assign the residue if possible
            if( isempty(RESIDUE_STRING) )
                RESIDUE_STRING = 'Resi';
            end
            obj.residue_array{Nc}= RESIDUE_STRING;
            obj.index_array{Nc}  = INDEX;
            obj.atom_array{Nc}   = ATOM_STRING;
            obj.curves{Nc}       = Curve;
            obj.Nc               = Nc;

            % Propagate properties to the new curve (and the rest)
            obj.updateCurves();
            
            if( OUTPUT_DEBUG )
                fprintf('\nVCPMG\tIntensity\teIntensity\tR2eff\teR2eff\n');
                disp([VCPMG' Intensity' eIntensity' R2eff' eR2eff']);                
                obj
            end
        end
        
        function allocateCurves( obj )
            %% Allocate a new instance of Curve for each curve in dataset
            % This should only be done ONCE, as soon as Nc are known
            % so subsequent references to curves are to single instance in memory
            for c = 1:obj.Nc
                obj.curves{c} = Curve;
            end
        end
        
       
        function calculateR2eff( obj )
            %% Calculate R2eff values from intensity and TCPMG
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_CALCULATE_R2EFF;            
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            if( ~isempty(obj.Intensity_matrix) && obj.TCPMG > 0 )
                                
                Nobs = size(obj.Intensity_matrix,2);                
                if( Nobs < length(obj.vcpmg) )
                    fprintf('\nAlter intensity matrix to remove first VCPMG=0 observations');
                    obj.Intensity_matrix
                    obj.eIntensity_matrix
                    
                    obj.Intensity_matrix = obj.Intensity_matrix(:,2:end);
                    obj.eIntensity_matrix = obj.eIntensity_matrix(:,2:end);
                    
                    obj.Intensity_matrix
                    obj.eIntensity_matrix                    
                end             
                
                obj.R2eff_matrix     = -log( obj.Intensity_matrix ) ./ obj.TCPMG;
                obj.eR2eff_matrix    = obj.eIntensity_matrix ./ ( obj.Intensity_matrix * obj.TCPMG );
                
                % Only keep the non-zero observations (vcpmg ~= 0)
                if( OUTPUT_DEBUG )                    
                    fprintf('\nRemoving vcpmg=0 elements');
                    fprintf('\nZero elements\n');
                    disp(find(obj.vcpmg == 0))
                    
                    fprintf('\nNon-zero elements\n');
                    disp(find(obj.vcpmg ~= 0))
                end
                
                oNonZeros               = find(obj.vcpmg ~= 0);
                obj.vcpmg               = obj.vcpmg(oNonZeros);
                
                obj.Intensity_matrix    = obj.Intensity_matrix(1:end, oNonZeros);
                obj.eIntensity_matrix   = obj.eIntensity_matrix(1:end, oNonZeros);
                
                obj.R2eff_matrix        = obj.R2eff_matrix(1:end, oNonZeros);
                obj.eR2eff_matrix       = obj.eR2eff_matrix(1:end, oNonZeros);
                obj.Nobs                = length(oNonZeros);
                
                if( OUTPUT_DEBUG )                    
                    fprintf('\nUpdated intensities and errors, after removing zero vcpmg points');
                    obj.vcpmg
                    obj.Intensity_matrix
                    obj.eIntensity_matrix
                    
                    fprintf('\nCalculated R2eff and errors, after removing zero vcpmg points');
                    obj.R2eff_matrix
                    obj.eR2eff_matrix
                    
                    fprintf('\nNumber of observations = %d', obj.Nobs);
                end
                
                obj.updateCurves();
                
            else
               error('R2eff requires defined Intensity and TCPMG'); 
            end            
        end
        
        function enforceMinimumError(obj, varargin)
            %% Enforce minimum error in data
            % Either supply minimum error (or use from parentSession)
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_ENFORCE_MIN_ERROR;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end
            
            if( nargin == 2 )
                MIN_F_ERROR = varargin{1};
            else
                MIN_F_ERROR = obj.parentSession.MIN_F_ERROR;
            end
            
            fError_matrix           = obj.eIntensity_matrix ./ obj.Intensity_matrix;
            iErrorTooSmall_matrix   = fError_matrix < MIN_F_ERROR;
            
                                    
            obj.eIntensity_matrix(iErrorTooSmall_matrix) = ...
                MIN_F_ERROR * obj.Intensity_matrix(iErrorTooSmall_matrix);
                        
            if( OUTPUT_DEBUG )
                %fError_matrix
                %iErrorTooSmall_matrix                
                if( ~isempty(iErrorTooSmall_matrix) )
                    fprintf('\n\tError too small (<%f) in %d observations', ...
                        MIN_F_ERROR, sum(sum(iErrorTooSmall_matrix)) );
                end
            end
            
            % Calculate R2eff if possible            
            obj.calculateR2eff()            
        end
        
        function exportData(obj, FILE)
            %% Export curves in dataset to a file
            % FILE = already opened file handle (1 => command window)
            
            
            fprintf(FILE, '\n%s', 'DATASET');
            fprintf(FILE, '\n%s\t%s', 'NAME',   obj.name);
            fprintf(FILE, '\n%s\t%s', 'AX',     obj.AX_String);
            fprintf(FILE, '\n%s\t%f', 'B0',     obj.B0);
            fprintf(FILE, '\n%s\t%f', 'TEMPC',  obj.Temp_C);
            fprintf(FILE, '\n%s\t%f', 'TCPMG',  obj.TCPMG);
            if( obj.SQX )
                SQX_String = 'TRUE';
            else
                SQX_String = 'FALSE';
            end
            fprintf(FILE, '\n%s\t%s', 'SQX',    SQX_String);
            fprintf(FILE, '\n%s', 'SETSPECS');
            
            %fprintf(FILE, '\n\n#%s\t%d', 'NumCurves',      obj.Nc);            
            
            % Export each curve
            for c = 1:obj.Nc                
                fprintf(FILE, '\n\n# Dataset %s, Curve %d/%d', obj.name, c, obj.Nc);
                
                curve = obj.curves{c};
                curve.exportData(FILE);                
            end
        end
                
        function specs_string = getSpecsString(obj)
            %% Return a string of specifications
            if( obj.SQX )
                SQ_string = 'SQ';
            else
                SQ_string = 'MQ';
            end
            specs_string = sprintf('%s-%dMHz-%dC', SQ_string, obj.B0, obj.Temp_C);            
        end
        
        function plot( obj, c )
            %%
            plot( obj.vcpmg, obj.R2eff_matrix(c,:), '-ok' )
            %title( sprintf('%s %d%s', obj.residue_array{c}, obj.index_array{c}, obj.atom_array{c}) )
            xlabel( '\nu_{CPMG} (Hz)' );
            ylabel( 'R_2^{eff} (Hz)' );
            
        end
        
        function readNlin( obj, nlin_filename, vcpmg_filename )
            % Read NMRPipe nlin.tab file, which contains RD curves for all NMR signals at a particular set of conditions
            % Ian Kleckner
            % 2008/05/06 Start code
            % 2008/05/05 Modified to read in a comment field in the second column
            % 2010/01/08 Change case of amino acid (e.g., TRP -> Trp)
            % 2011/01/12 Change "data" -> "dataset"
            % 2011/01/18 Update for classes
            % 2011/06/05 Update for arbitrary lines before FORMAT token
            %
            % INPUT
            %  nlin.tab data file name
            %  vcpmgfile with values of vcpmg (vCPMG or time,etc.) delimited by whitespaces
            %  (Optional) sequence file for AA sequence
            %
            % OUTPUT
            %  dataset (from Dataset class)
            %
            % TO DO
            %
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_READ_NLIN_FILE;
            if( OUTPUT_DEBUG )
                [ST,I] = dbstack;                
                fprintf('\n\nFUNCTION: %s', ST(1).name);
            end

            % ================= Read the nlin.tab file =========================== %            
            if( ~exist(nlin_filename, 'file') )
                error('Cannot find nlin file %s', nlin_filename);
            end
            if( ~exist(vcpmg_filename, 'file') )
                error('Cannot find vcpmg file %s', vcpmg_filename);
            end
            
            FILE    = fopen(nlin_filename, 'r');
            lines   = textscan(FILE, '%s', 'delimiter', '\n');
            Nlines  = length(lines{1});
            fclose(FILE);

            % TODO % Change this code for more robust reading
            
            % Read input lines until finding FORMAT token
            for l = 1:Nlines
                % Get the line string
                line = lines{1}(l);
                % Convert line from CELL to CHAR
                while( iscell(line) )
                    line = line{1};
                end
                
                % Skip the first lines of comments
                %line = textscan(FILE, '%s',1,'delimiter','\n');
                                                
                if( isempty(line) )
                    Ntokens     = 0;
                    firstToken  = '';
                    
                else
                    tokens  = textscan(char(line), '%s');
                    Ntokens = size(tokens{1},1);
                    
                    firstToken = tokens{1}(1);
                    firstToken = firstToken{1};
                end
                
                if( strcmpi(firstToken,'FORMAT') )
                    % Found "TOKEN"
                    if( OUTPUT_DEBUG )
                        fprintf('\n\tFound "TOKEN" line at line number %d', l);
                    end
                    break
                end
            end
            currentLine = l;
            
            if( currentLine == Nlines )
                error('FORMAT string not found');
            end
            
            % The index at which each key variable is found
            % TODO % NLIN.tab files may differ here in different versions of NMRPipe
            token_index = 1;
            token_atom  = 23;
            token_data  = 26;            
                        
            % Now "line" holds the line for the FORMAT information

            % Count the number of spaces in this line of text
            %numSpaces = size(find(char(inputText{1}(1)) == ' '),2);
            numSpaces = size(find(line == ' '),2);

            % There are token_data number of fields which are not series data poinst Z_A0 Z_A1, ... Z_An
            Nobs = numSpaces+1 - token_data;
            
            if( OUTPUT_DEBUG )
                fprintf('\n\tFound %d observations in NLIN file', Nobs);
                fprintf('\n\tReading data points...');
            end
            
            
            % Read input lines until finding a line with DATA
            for l = currentLine+1:Nlines
                % Get the line string
                line = lines{1}(l);
                % Convert line from CELL to CHAR
                while( iscell(line) )
                    line = line{1};
                end
                
                % Skip the first lines of comments
                %line = textscan(FILE, '%s',1,'delimiter','\n');
                                                
                if( isempty(line) )
                    Ntokens     = 0;
                    firstToken  = '';
                    
                else
                    tokens  = textscan(char(line), '%s');
                    Ntokens = size(tokens{1},1);
                end
                
                if( Ntokens > Nobs )
                    % Found a data line
                    if( OUTPUT_DEBUG )
                        fprintf('\n\tFound a data line at line number %d', l);
                    end
                    break
                end
            end
            currentLine = l;
            

            % Read through each line of the file (to get all peaks)
            peak=1;
            for l = currentLine:Nlines                
                % Get the line string
                line = lines{1}(l);
                % Convert line from CELL to CHAR
                while( iscell(line) )
                    line = line{1};
                end
                
                if( OUTPUT_DEBUG )
                    fprintf('\n\tReading line %d/%d: "%s"...', l, Nlines, line);
                end
                
                if( ~isempty(line) )                                    
                    % Read the first set of data for the peak (shown below in table)                    
                    tokens = textscan(line, '%s');
                    tokens = tokens{1};

                    if( OUTPUT_DEBUG )
                        fprintf('\n\t\tRead text:');
                    end

                    % If there is no more data (i.e. EOF) then exit the loop
                    if( ~isempty(tokens) )                        
                        % Read index and atom
                        obj.index_array{peak}   = str2double( tokens{token_index} );
                        obj.atom_array{peak}    = tokens{token_atom};

                        % Read the Z_A0, Z_A1, ... Z_A(Nobs-1) values
                        for o = 1:Nobs
                            %inputText = textscan(line, '%f');
                            Intensity_Matrix(peak,o) = str2double( tokens{token_data+o-1} );
                        end
                        
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tRead index=%d, atom=%s, and Intensity_Matrix=(%d points) for peak=%d:', ...
                                obj.index_array{peak}, obj.atom_array{peak}, Nobs, peak);
                        end

                        peak = peak+1;
                       
                    else
                        if( OUTPUT_DEBUG )
                            fprintf('\n\t\tNo data to parse in input line:');
                            inputText
                        end
                    end
                   
                else
                    if( OUTPUT_DEBUG )
                        fprintf('\n\t\tEmpty');
                    end
                end
            end
            numPeaks = peak-1;
            
            %% ================= Read the vcpmg file =========================== %
            % Now read in the vcpmg file which labels each of the Z_A0, Z_A1, ...
            % This can be with Vcpmg or with time, etc.
            FILE = fopen(vcpmg_filename, 'r');
            inputText = textscan(FILE, '%f');

            VCPMG = inputText{1};
            fclose(FILE);
            
            %% ================= Store data to class =========================== %
            obj.nlin_filename   = nlin_filename;
            obj.vcpmg_filename  = vcpmg_filename;
            obj.Nc              = numPeaks;            
            obj.allocateCurves();
            
            if( strcmp(obj.name,'') )
                obj.name = nlin_filename;
            end
            
            if( isempty(obj.residue_array) )
                for c = 1:obj.Nc
                    obj.residue_array{c} = 'Resi';
                end
            end
            
            if( OUTPUT_DEBUG )
                fprintf('\nBefore removing duplicate entries');
                VCPMG
                Intensity_Matrix
            end

            % Calculate the errors in Intensity values
            [ obj.vcpmg, obj.Intensity_matrix, obj.eIntensity_matrix] = ...
                Dataset.calculateErrorsUsingDuplicates( VCPMG, Intensity_Matrix );
            
            if( OUTPUT_DEBUG )
                fprintf('\nAfter calculating errors, removing duplicates');
                obj.vcpmg
                obj.Intensity_matrix
                obj.eIntensity_matrix
            end
            
            % Enforce minimum error in data
            obj.enforceMinimumError();
            
            % Calculate R2eff if possible
            if( obj.TCPMG > 0 )
                obj.calculateR2eff()
            end
        end
            
              
        function readR2Table( obj, R2table_filename )
            %%
            
            %obj.curves{c} = Curve;
        end
        
        function removeCurve( obj, curve )
           %% Remove this curve from the dataset (and all other instances)
           % TODO % Finish removeCurve()
            
        end
        
        function saveDataFromTable( obj, table_data )
            %% Set parameters using values from table
            
            % Save each property, row by row
            %  Column 1 contains property name
            %  Column 2 contains value entered into the table
            for row = 1:size(table_data,1)
                propety_name = table_data{row,1};
                eval( sprintf('obj.%s = table_data{row,2};', propety_name) );                
            end
            % Convert Celcius to Kelvin
            obj.Temp = obj.Temp_C + 273;
            
            obj.updateCurves()
        end
        
        function setParentSession(obj, parentSession)
            %% Set the parent session for this dataset
            obj.parentSession = parentSession;
        end
        
        function setSequence( obj, sequence_array )
            %% Set residue names
            OUTPUT_DEBUG = obj.parentSession.OUTPUT_DEBUG_SET_SEQUENCE;
            
            for c = 1:obj.Nc
                if( obj.index_array{c} <= length(sequence_array) )
                    obj.residue_array{c} = sequence_array{obj.index_array{c}};
                else
                    obj.residue_array{c} = 'Resi';
                end
                
                if( OUTPUT_DEBUG )
                    fprintf('\n\tSet curve #%d at index %d with residue %s', ...
                        c, obj.index_array{c}, obj.residue_array{c});
                end
            end
            
            % Propagate this to the curves
            obj.updateCurves();
        end
        
        function obj = setSpecs( obj, name, AX_String, B0, Temp, TCPMG, SQX )            
            %% Set the specifications for the dataset (name=[] => autoname)
            % See text above for descriptions of properties            
            obj.AX_String   = AX_String;
            obj.B0          = B0;
            obj.Temp        = Temp;
            obj.Temp_C      = Temp-273;
            obj.TCPMG       = TCPMG;
            obj.SQX         = SQX;
            
            % Autoname if desired
            if( isempty(name) )
                % Automatically set name (e.g., 800-MQ-25C)
                if( obj.SQX )
                    QC_String = 'SQ';
                else
                    QC_String = 'MQ';
                end
                obj.name = sprintf('%d-%s-%dC', round(obj.B0), QC_String, obj.Temp-273);
            else
                obj.name    = name;
            end
            
            % Calculate R2eff if possible
            if( ~isempty(obj.Intensity_matrix) )
                obj.calculateR2eff();
                
            else                
                obj.updateCurves()
            end
        end
        
        function sortCurves( obj )
            %% Sort the curve_array by index (?)
            % TODO % Must also sort other _arrays that index by "c"
            % TODO % Perhaps this should go to higher-order class / function?
            % TODO % Also may want to sort curves in curveset, curveset in group, etc.
            
            % Sort the indices in the groups
            cIndex = zeros(1,obj.Nc);
            for c = 1:obj.Nc
               cIndex(c) = obj.curves{c}.index; 
            end            
            [VOID, c_sort_array] = sort(cIndex);
            
            % Use that sorted array to rearrange order of groups
            new_curves = cell(1,obj.Nc);
            for c_sort_index = 1:obj.Nc
                % The index of csort which should go next
                new_curves{c_sort_index} = obj.curves{ c_sort_array(c_sort_index) };
            end
            
            obj.curves = new_curves;
        end
                    
        function updateCurves( obj )
            %% Update parameters for each RD curve in the dataset
            for c = 1:obj.Nc               
                if( ~isempty(obj.index_array) )
                    obj.curves{c}.setAssignment( obj.index_array{c}, obj.atom_array{c}, obj.residue_array{c} );
                end
                if( obj.B0 > 0 )
                    obj.curves{c}.setConditions( obj.B0, obj.Temp, obj.TCPMG, obj.SQX, obj.AX_String );
                end
                if( ~isempty(obj.vcpmg) )
                    obj.curves{c}.setData( obj.vcpmg, obj.R2eff_matrix(c,:), obj.eR2eff_matrix(c,:) );
                end
                
                obj.curves{c}.setDataset( obj );
            end
        end
    end
    
    methods (Static = true)
        function [newVCPMG, newIntensity_Matrix, eIntensity_Matrix] = ...
                calculateErrorsUsingDuplicates( VCPMG, Intensity_Matrix )
            %% Calculate the errors in Intensity values using duplicate measures
            % VCPMG is a 1D array of X-values (column vector)
            % Intensity_Matrix is an nD array of Y-values, each column a different dataset
            % Check if there are any duplicate vcpmg values recorded
            
            % Ensure VCPMG is a ROW vector
            %  100 200 300 400 400 500 ...
            [Nrows, Ncols] = size(VCPMG);
            if( Ncols == 1 )
                VCPMG = VCPMG';
                [Nrows, Ncols] = size(VCPMG);
            end
            Nobs = Ncols;
            
            % Ensure Intensity_Matrix is a ROW vector (each column a VCPMG value)
            %  I1     I2      I3 ...
            %  I1     I2      I3 ...
            [Nrows, Ncols] = size(Intensity_Matrix);
            
            % Only one column may mean the matrix is transposed => fix it
            if( Ncols == 1 )
                Intensity_Matrix = Intensity_Matrix';
                [Nrows, Ncols] = size(Intensity_Matrix);
            end
            
            if( Ncols ~= Nobs )
                error('Must provde same number of Intensities as VCPMG values');
            end
            numPeaks = Nrows;

            % For each value of vcpmg in the VCPMG array, search the VCPMG array for the value
            dupVcpmgIndex   = 1;
            
            % Holds indices of each set of duplicate values
            dupVcpmgs       = [];
            
            for i = 1:Nobs
                % Check if there are other occurrances of the value vcpmg
                vcpmg = VCPMG(i);
                dupVcpmg = find( VCPMG == vcpmg );
                
                % If there is a duplicate value, then note the pair
                if( size(dupVcpmg,2) > 1 )
                    dupVcpmgs(dupVcpmgIndex,:) = dupVcpmg;
                    dupVcpmgIndex = dupVcpmgIndex + 1;
                end
            end

            % Now get only the unique values
            dupVcpmgs = unique(dupVcpmgs,'rows');

            % Proceed only if there are duplicate values
            if( size(dupVcpmgs,2) > 1 )

                % Calculate errors for each peak
                for peak = 1:numPeaks
                    uIntensity_Matrix  = [];
                    sIntensity_Matrix  = [];
                    feIntensity_Matrix = [];

                    % For each set of duplicate values
                    for j = 1:size(dupVcpmgs,1)

                        % Find the duplicate entries
                        dupVcpmg = dupVcpmgs(j,:);

                        % Calcualte the mean
                        sum = 0;
                        for i = 1:size(dupVcpmg,2)
                            sum = sum + Intensity_Matrix(peak,dupVcpmg(i));
                        end
                        uIntensity_Matrix = sum/i;

                        % And the standard deviation
                        ss = 0;
                        for i = 1:size(dupVcpmg,2)
                            ss = ss + (uIntensity_Matrix - Intensity_Matrix(peak,dupVcpmg(i)))^2;
                        end
                        sIntensity_Matrix = sqrt( ss/(i-1) );

                        % Calcualte the fractional error
                        feIntensity_Matrix(j) = sIntensity_Matrix / uIntensity_Matrix;

                        % Set the value of EACH duplicate to its mean
                        for i = 1:size(dupVcpmg,2)
                            newIntensity_Matrix(peak,dupVcpmg(i))   = uIntensity_Matrix;
                            newVCPMG(dupVcpmg(i))                   = VCPMG(dupVcpmg(i));
                        end
                    end

                    % Copy the rest of the values over to newIntensity_Matrix, skipping other
                    % duplicate entries
                    for j = 1:Nobs
                        if( isempty( find(dupVcpmgs == j) ) )
                            newIntensity_Matrix(peak,j) = Intensity_Matrix(peak,j);
                            newVCPMG(j) = VCPMG(j);
                        end
                    end

                    % The fractional error for all intensities is the average
                    % of all the fractional errors for each group of duplicates
                    fError = mean(feIntensity_Matrix(:));

                    % Calculate the errors in Intensity_Matrix using fractional error
                    eIntensity_Matrix(peak,:) = newIntensity_Matrix(peak,:).*fError;
                end

                % Now delete ALL BUT ONE of each of the duplicate columns from
                % VCPMG, newIntensity_Matrix and eIntensity_Matrix
                % Find the positions to delete
                idelete = 0;
                for i = 1:size(dupVcpmgs,1)
                    % Iterating from 2...size(dupVcpmgs,2) preserves the first repeat value
                    for j = 2: size(dupVcpmgs,2)
                        idelete = idelete+1;
                        deleteMe(idelete) = dupVcpmgs(i,j);
                    end
                end

               % Delete them from largest index to smallest index
               % This way, subsequent deletion indices remain accurate
               deleteMe = sort( deleteMe, 'descend' );
               for i = 1:idelete
                    newVCPMG(deleteMe(i))               = [];
                    newIntensity_Matrix(:, deleteMe(i)) = [];
                    eIntensity_Matrix(:, deleteMe(i))   = [];
               end
               
            else                
                % In the absence of duplicates, errors are unknown
                newIntensity_Matrix     = Intensity_Matrix;
                eIntensity_Matrix       = NaN .* Intensity_Matrix;
                newVCPMG                = VCPMG;
            end
        end
    end
end

