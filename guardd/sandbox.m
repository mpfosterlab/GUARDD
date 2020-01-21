% File for testing assorted GUARDD functions
%
% (C) Ian Kleckner [ian.kleckner@gmail.com]
%  Foster Lab, The Ohio State University
% GUARDD software [http://code.google.com/p/guardd/]
%  GNU GPL3 License
%

clear all
clc

%% 2011/09/02 Perform analyzeMe() on all fitResult instances
% Re-compute the resultsMatrix for PhiexX
for g = 1:session.Ng
    fprintf('\nWorking on group %d/%d', g, session.Ng);
    group = session.groups{g};
    for f = 1:group.Nf
        fitResult = group.fitResults{f};        
        fitResult.analyzeMe();
    end
    
    for fG = 1:length(group.fitResults_Grid)
        fprintf('\n\tWorking on grid %d/%d', fG, length(group.fitResults_Grid));
        fitResult = group.fitResults_Grid{fG};
        fitResult.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_NoEx) )
        group.fitResult_NoEx.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_Best) )
        group.fitResult_Best.analyzeMe();    
    end
end

%% Testing subplots (2011/06/14)

TOP_SIZE = 4;

h = subplot(TOP_SIZE,1,1:TOP_SIZE-1);
hold(h, 'all');

X = 1:1:100;
plot(h, X, sin(X));
plot(h, X, cos(X));


h = subplot(TOP_SIZE,1,TOP_SIZE);
plot(h, X, sin(X));

% Can I go back to subplot 1 without resettting it? => Yes!
h = subplot(TOP_SIZE,1,1:TOP_SIZE-1);
plot(h, X, 2*sin(X));

%% Get symbol and color (2011/06/10)
% FOR: function [symbolChar, colorRGB] = getPlotSymbolAndColor(plotNum)

plotNum = 1:10;
Ncolors = 3;



%% 2011 / 06 / 01 - PErform rate analysis on all fits
% Convert kcal -> cal
for g = 1:session.Ng
    fprintf('\nWorking on group %d/%d', g, session.Ng);
    group = session.groups{g};
    for f = 1:group.Nf
        fitResult = group.fitResults{f};        
        fitResult.rateAnalysis.analyzeMe();
    end
    
    for fG = 1:length(group.fitResults_Grid)
        fprintf('\n\tWorking on grid %d/%d', fG, length(group.fitResults_Grid));
        fitResult = group.fitResults_Grid{fG};
        fitResult.rateAnalysis.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_NoEx) )
        group.fitResult_NoEx.rateAnalysis.analyzeMe();
    end
    
    if( ~isempty(group.fitResult_Best) )
        group.fitResult_Best.rateAnalysis.analyzeMe();    
    end
end

%% 2011 / 06 / 05 - Reset parameter display for dwHppm and dwXppm
% to be shown in the group display
session.resetParamDisplay();


%% 2011 / 06 / 01 - PErform rate analysis on all best fits
for g = 1:session.Ng
    group = session.groups{g};
    group.fitResult_Best.rateAnalysis.analyzeMe();
end

%% 2011 /05 / 31 - Check use of "continue" command
%{
for c = 1:10
    fprintf('\n\nc = %d',c);
    fprintf('\n\tPre-check');
    if( c == 3 )
        fprintf('\n\tAborting loop...');
        continue;
        fprintf('\n\tPost-continue');
    end
    fprintf('\n\tPassed');
end
%}


%% 2011 / 05 / 26 - Enforce minimum errors in data
%{
for ds = 1:session.Nds
    fprintf('\nWorking on dataset %d, %s', ds, session.datasets{ds}.name);
    session.datasets{ds}.enforceMinimumError(session.MIN_F_ERROR);
end
%}


%% 2011/05/22 - Convert group notes
%{
session.groups{get(handles.listbox_g,'Value')}

for g = 1:session.Ng
    group = session.groups{g};
    
    group.setNote( group.notes );
end
%}


%% Correct stereo-assignments (2011 / 05 / 22)
% Put this in GUARDD.m -> check_variables()
%{
OLD_LEU_DELTA_1     = '\delta_1';
NEW_LEU_DELTA_1     = '\delta_2';

OLD_LEU_DELTA_2     = '\delta_2';
NEW_LEU_DELTA_2     = '\delta_1';

OLD_VAL_GAMMA_1     = '\gamma_1';
NEW_VAL_GAMMA_1     = '\gamma_2';

OLD_VAL_GAMMA_2     = '\gamma_2';
NEW_VAL_GAMMA_2     = '\gamma_1';

fprintf('\n\nWORKING ON DATASETS');

% Rename all curves in datasets
for ds = 1:session.Nds
    dataset = session.datasets{ds};
    fprintf('\n\n(%d / %d) Working on %s', ds, session.Nds, dataset.name);
    
    for c = 1:dataset.Nc
        curve = dataset.curves{c};
        old_atom = curve.atom;
        
        switch upper( curve.residue )
            case 'LEU'
                switch( old_atom )
                    case OLD_LEU_DELTA_1
                        new_atom = NEW_LEU_DELTA_1;
                    case OLD_LEU_DELTA_2
                        new_atom = NEW_LEU_DELTA_2;
                    otherwise
                        new_atom = '???';
                end
            case 'VAL'
                switch( old_atom )
                    case OLD_VAL_GAMMA_1
                        new_atom = NEW_VAL_GAMMA_1;
                    case OLD_VAL_GAMMA_2
                        new_atom = NEW_VAL_GAMMA_2;
                    otherwise
                        new_atom = '???';
                end
                
            case 'ILE'
                % Nothing to do
                new_atom = curve.atom;
                    
            otherwise
        end
        
        fprintf('\n\t(%d / %d) Working on %s', c, dataset.Nc, curve.name);
        curve.setAssignment( curve.index, new_atom, curve.residue )        
        fprintf(' -> %s', curve.name);
    end
end


fprintf('\n\nWORKING ON GROUPS');
for g = 1:session.Ng
    group = session.groups{g};
    
    % Find the residue via name
    group_residue = group.name(1:3);
    old_name = group.name;
    
    switch upper( group_residue )
            case 'LEU'
                if( ~isempty(strfind(old_name, OLD_LEU_DELTA_1)) )
                    old_atom = OLD_LEU_DELTA_1;
                    new_atom = NEW_LEU_DELTA_1;
                    
                elseif( ~isempty(strfind(old_name, OLD_LEU_DELTA_2)) )
                    old_atom = OLD_LEU_DELTA_2;
                    new_atom = NEW_LEU_DELTA_2;
                else
                    '??';
                end
                
            case 'VAL'
                if( ~isempty(strfind(old_name, OLD_VAL_GAMMA_1)) )
                    old_atom = OLD_VAL_GAMMA_1;
                    new_atom = NEW_VAL_GAMMA_1;
                    
                elseif( ~isempty(strfind(old_name, OLD_VAL_GAMMA_2)) )
                    old_atom = OLD_VAL_GAMMA_2;
                    new_atom = NEW_VAL_GAMMA_2;
                else
                    '??';
                end
                
            case 'ILE'
                % Nothing to do
                old_atom = 'X';
                new_atom = old_atom;
                    
            otherwise
    end
    
    fprintf('\n\n(%d / %d) Working on %s', g, session.Ng, group.name);
   
    % Sometimes the name is formatted differently than default
    new_name = strrep(old_name, old_atom, new_atom);
    group.setName(new_name);
    fprintf(' -> %s', group.name);
        
    
    for cs = 1:group.Ncs
        curveset = group.curvesets{cs};        
        old_name = curveset.name;
        old_atom = curveset.atom;
        
        switch upper( curveset.residue )
            case 'LEU'
                switch( old_atom )
                    case OLD_LEU_DELTA_1
                        new_atom = NEW_LEU_DELTA_1;
                        
                    case OLD_LEU_DELTA_2
                        new_atom = NEW_LEU_DELTA_2;
                        
                    otherwise
                        new_atom = '???';
                end
            case 'VAL'
                switch( old_atom )
                    case OLD_VAL_GAMMA_1
                        new_atom = NEW_VAL_GAMMA_1;
                        
                    case OLD_VAL_GAMMA_2
                        new_atom = NEW_VAL_GAMMA_2;
                        
                    otherwise
                        new_atom = '???';
                end
                
            case 'ILE'
                % Nothing to do
                new_atom = curveset.atom;
                    
            otherwise
        end
        
        fprintf('\n\t(%d / %d) Working on %s', cs, group.Ncs, curveset.name);
        
        % Do NOT change curves within, because those have already been
        % changed above
        curveset.setAtom_DoNotChangeCurves( new_atom );
        
        % Sometimes the name is formatted differently than default
        new_name = strrep(old_name, old_atom, new_atom);
        curveset.setName(new_name);
        fprintf(' -> %s', curveset.name);
    end
end
%}

%% 2011/05/20
%{
X = rand(50,1);
Y = rand(10,1);

%f = figure;
f = 1;
clf(f);

h = axes;
cla(h);

hold(h,'all');

[n,xout] = hist(h,X);
bar(xout,n, 'FaceColor', 'g', 'EdgeColor', 'k');
%hh = findobj(h,'Type','patch');
%set(hh, 'FaceColor', 'g', 'EdgeColor', 'k');

[n,xout] = hist(h,Y);
bar(xout,n, 'FaceColor', 'r', 'EdgeColor', 'k');
%}


%% 2011/05/19
% 2011/06/09 -> This is a deprecated use of addData (now requires INDEX, ATOM, RESIDUE, ...)
%{
d = Dataset;
d.setSpecs('NAME', '15N', 800, 298, 0.02, true)

VCPMG = [100 200 300 400 500];
R2 = [50 45 40 35 35];
ER2 = [5 4 3 3 2];

INDEX = 10;
ATOM = 'NH';

d.addData(INDEX, ATOM, VCPMG, R2, ER2, 'R2')

d.addData(INDEX+1, ATOM, VCPMG, R2.*1.1, ER2, 'R2')

VCPMG = [0 100 200 300 400 500];
I = [1 .5 .6 .7 .8 .8];
EI = [0 .1 .1 .1 .1 .1];

d.addData(INDEX+2, ATOM, VCPMG, I, EI, 'INTENSITY')


VCPMG = [0 100 200 300 400 500 500];
I = [1 .5 .6 .7 .8 .8 .9];
EI = [];

d.addData(INDEX+3, ATOM, VCPMG, I, EI, 'INTENSITY')
%}

