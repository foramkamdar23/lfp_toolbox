%% FEELpos vs FEELneu
disp('--- pos vs neu ---');

if isfield(stat_freq_pos_vs_neu, 'posclusters') && ~isempty(stat_freq_pos_vs_neu.posclusters)
    p_pos_valpos = stat_freq_pos_vs_neu.posclusters(1).prob;
    fprintf('Largest POSITIVE cluster p = %.4f\n', p_pos_valpos);
else
    fprintf('No positive clusters found.\n');
end

if isfield(stat_freq_pos_vs_neu, 'negclusters') && ~isempty(stat_freq_pos_vs_neu.negclusters)
    p_neg_valpos = stat_freq_pos_vs_neu.negclusters(1).prob;
    fprintf('Largest NEGATIVE cluster p = %.4f\n', p_neg_valpos);
else
    fprintf('No negative clusters found.\n');
end

%% FEELneg vs FEELneu
disp('--- neg vs neu ---');

if isfield(stat_freq_neg_vs_neu, 'posclusters') && ~isempty(stat_freq_neg_vs_neu.posclusters)
    p_pos_valneg = stat_freq_neg_vs_neu.posclusters(1).prob;
    fprintf('Largest POSITIVE cluster p = %.4f\n', p_pos_valneg);
else
    fprintf('No positive clusters found.\n');
end

if isfield(stat_freq_neg_vs_neu, 'negclusters') && ~isempty(stat_freq_neg_vs_neu.negclusters)
    p_neg_valneg = stat_freq_neg_vs_neu.negclusters(1).prob;
    fprintf('Largest NEGATIVE cluster p = %.4f\n', p_neg_valneg);
else
    fprintf('No negative clusters found.\n');
end