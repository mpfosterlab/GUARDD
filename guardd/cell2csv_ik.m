function cell2csv_ik(datName,cellArray,seperator)
% Writes cell array content into a *.csv file.
%
% CELL2CSV(datName,cellArray,seperator,excelVersion)
% Modified by Ian Kleckner on 03/02/2010
% 2011/06/10 Modify file to append to current contents
%
% datName      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% seperator    = seperating sign, normally:',' (it's default)
% excelVersion = depending on the Excel Version, the cells are put into
%                quotes before added to the file (only numeric values)
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler
%
% 2010/03/02 IK: Modified by Ian Kleckner (IK)
% 2011/06/05 IK: Replace "%" with "%%" for proper export

if seperator ~= ''
    seperator = ',';
end

% Removed by IK
%{
if excelVersion > 2000
    seperator = ';';
end
%}

% Debugging output?
output_debug = false;

% Append the file
FILE = fopen(datName,'a');

for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)
        
        var = eval(['cellArray{z,s}']);
        
        %{
        if size(var,1) == 0
            var = '';
        end
        %}
        
        % VVV "if" Added by IK on 03/02/2010
        if( size(var,1) > 0 )
        
            if isnumeric(var) == 1
                var = num2str(var);
            end

            if islogical(var) == 1
                if var == 1
                    var = 'TRUE';
                else
                    var = 'FALSE';
                end
            end

            % Added by Ian Kleckner 03/02/2010
            if( iscell(var) )
                var = var{1};
            end

            % Do this for *every* var (IK)
            var = ['"' var '"'];
            
            % Replace '\' char with '-'
            var = strrep(var,'\','-');
            
            % Replace % with %% (ik, 2011/06/05)
            var = strrep(var,'%','%%');

            fprintf(FILE,var);
            
            
            if( output_debug )
                fprintf('%s\t',var);
            end

            % If its not the last string to output, add the seperator
            if s ~= size(cellArray,2)
                fprintf(FILE,seperator);
            end
        end
        % ^^^ "if" Added by IK on 03/02/2010
    end
    fprintf(FILE,'\n');
    
    if( output_debug )
        fprintf('\n');
    end
end
fclose(FILE);

if( output_debug )
    fprintf('\n\nIn cell2csv_ik\tRows = %d, Cols = %d', size(cellArray,1), size(cellArray,2) )
end

%{
UNMODIFIED VERSION

function cell2csv(datName,cellArray,seperator,excelVersion)
% Writes cell array content into a *.csv file.
% 
% CELL2CSV(datName,cellArray,seperator,excelVersion)
%
% datName      = Name of the file to save. [ i.e. 'text.csv' ]
% cellarray    = Name of the Cell Array where the data is in
% seperator    = seperating sign, normally:',' (it's default)
% excelVersion = depending on the Excel Version, the cells are put into
%                quotes before added to the file (only numeric values)
%
%         by Sylvain Fiedler, KA, 2004
% updated by Sylvain Fiedler, Metz, 06
% fixed the logical-bug, Kaiserslautern, 06/2008, S.Fiedler

if seperator ~= ''
    seperator = ',';
end

if excelVersion > 2000
    seperator = ';';
end

datei = fopen(datName,'w');

for z=1:size(cellArray,1)
    for s=1:size(cellArray,2)
        
        var = eval(['cellArray{z,s}']);
        
        if size(var,1) == 0
            var = '';
        end
        
        if isnumeric(var) == 1
            var = num2str(var);
        end
        
        if islogical(var) == 1
            if var == 1
                var = 'TRUE';
            else
                var = 'FALSE';
            end
        end
        
        if excelVersion > 2000
            var = ['"' var '"'];
        end
        fprintf(datei,var);
        
        if s ~= size(cellArray,2)
            fprintf(datei,seperator);
        end
    end
    fprintf(datei,'\n');
end
fclose(datei);

%}
