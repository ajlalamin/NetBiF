%% NetBiF - Network Bifurcation Finder
% AL JAY LAN ALAMIN (ajalamin)
% Date created: 23 July 2025
% Last Modified: 26 August 2026

% clc;
rng(2);   % Reproducible randomisation

inspect_model(model)
[bif_species, bif_rxns, bif_sources, ...
        lbs_list, model, amat, A_sym, smat, kmat] = bifurcation_finder(model);

%% ===================================================================
%  MAIN FUNCTION: bifurcation_finder
% ===================================================================
function [bif_species, bif_rxns, bif_sources, ...
        lbs_list, model, amat, A_sym, smat, kmat] = bifurcation_finder(model)

    fprintf('\n========================================\n');
    fprintf('  NetBiF: Network Bifurcation Finder\n');
    fprintf('========================================\n\n');
    fprintf('Model: %s\n\n', model.id);

    tic;

    % == STEP 0: Preprocessing ==========================================
    fprintf('Step 0: Preprocessing model...\n');
    [model, m] = modelSpecies(model);
    [smat, kmat, rmat, ~, r] = stoichMatrix(model, m);
    fprintf('  Network: %d species, %d reactions\n\n', m, r);

    % == STEP 1: Augmented matrix A =====================================
    fprintf('Step 1: Computing augmented matrix A...\n');
    [amat, ~, ~] = compute_amat(smat, kmat, rmat, r, m);

    if abs(det(amat)) < eps
        fprintf('  WARNING: A is singular (det A = %g)\n', det(amat));
        ainv = pinv(amat);
    else
        ainv = inv(amat);
    end
    fprintf('  Matrix dimensions: (%d × %d)\n\n', size(amat, 1), size(amat, 2));

    % == STEP 2: Labeled Buffering Structures ===========================
    fprintf('Step 2: Identifying Labeled Buffering Structures (LBS)...\n');
    lbs_list = enumerate_lbs(amat, ainv, r, m);
    display_lbs(lbs_list, model);

    % == STEP 3: Symbolic A =============================================
    fprintf('\nStep 3: Creating symbolic representation...\n');
    A_sym = symbolify_amat(amat, model, r, m);

    % == STEP 4: Determinant structures =================================
    fprintf('\nStep 4: Building determinant structures (Γβ)...\n');
    det_struct_list = build_determinant_structures(lbs_list, model, m, r);

    % == STEP 5: Detect bifurcations ====================================
    fprintf('\nStep 5: Analyzing bifurcation candidates...\n');
    [bif_species, bif_rxns, bif_sources] = analyze_bifurcations( ...
        det_struct_list, A_sym, model);

    % == STEP 6: Summary ================================================
    display_bifurcation_summary(bif_species, bif_rxns, bif_sources, model, lbs_list);

    toc;

    NetBiF_classify(bif_species, bif_rxns, bif_sources, ...
        lbs_list, model, amat, A_sym, smat, kmat);
end


%% ===================================================================
%  CORE ALGORITHMS
% ===================================================================

function [amat, cmat, dmat] = compute_amat(smat, kmat, rmat, r, m)
% COMPUTE_AMAT - Build the augmented Jacobian matrix A.
%
% Layout:
%   A = [ drdx  |  C' ]    (r × m)  |  (r × c)
%       [ -D    |  0  ]    (d × m)  |  (d × c)
%
% drdx: randomised partial derivatives ∂rᵢ/∂xⱼ.
%   Nonzero entries are positions where species j is a substrate of
%   reaction i (kmat(j,i) > 0 or rmat(j,i) > 0). 
%
% C = null(smat,  'rational')   - kernel of stoichiometric matrix (cycles)
% D = null(smat', 'rational')   - cokernel (conservation laws)

    drdx = randomize_nonzero_entries((rmat + kmat)' + rmat');  % r × m

    cmat = null(smat,  'rational')';   % c × r
    dmat = null(smat', 'rational')';   % d × m

    c = size(cmat, 1);
    d = size(dmat, 1);

    % assert(m + c == r + d, 'Fredholm theorem violated: m+c ≠ r+d');

    amat = zeros(r + d, m + c);
    amat(1:r,     1:m)     =  drdx;    % top-left  : ∂r/∂x
    amat(1:r,     m+1:end) =  cmat';   % top-right : C'
    amat(r+1:end, 1:m)     = -dmat;    % bot-left  : -D
end


function lbs_list = enumerate_lbs(amat, ainv, r, m)
% ENUMERATE_LBS - Find all Labeled Buffering Structures.
%
% Method: expand candidate subnetworks from each reaction, filter for
% output-completeness, then keep those with index = |m| - |n| + nullity = 0.

    lbs_list = {};

    smat_resp = ainv(1:m, :);
    smat_bin  = abs(smat_resp) > eps;

    % Build candidate subnetworks
    affected_subnetworks = cell(1, r);
    for i = 1:r
        affected_subnetworks{i} = expand_bs_subnetwork(smat_bin, i);
    end

    n_subnets        = numel(affected_subnetworks);
    species_matrix   = false(m, n_subnets);
    reactions_matrix = false(r, n_subnets);

    for i = 1:n_subnets
        species_matrix(affected_subnetworks{i}.xs,  i) = true;
        reactions_matrix(affected_subnetworks{i}.rs, i) = true;
    end

    % Filter: output-complete (no regulated reaction left outside)
    rmat_top = amat(1:r, 1:m);
    output_complete = true(n_subnets, 1);
    for a = 1:n_subnets
        xs = find(species_matrix(:, a));
        rs = find(reactions_matrix(:, a));
        if isempty(xs) || isempty(rs), output_complete(a) = false; continue; end
        regulated_r = any(rmat_top(:, xs) ~= 0, 2);
        missing = setdiff(find(regulated_r), rs);
        if ~isempty(missing), output_complete(a) = false; end
    end

    if sum(~output_complete) > 0
        fprintf('  Filtered %d subnetworks (not output-complete)\n', sum(~output_complete));
    end

    % Check index = 0 (buffering structure condition)
    for i = 1:n_subnets
        if ~output_complete(i), continue; end

        xs_real = affected_subnetworks{i}.xs(:)';
        rs_real = affected_subnetworks{i}.rs(:)';
        xs_real = xs_real(xs_real <= m);
        rs_real = rs_real(rs_real <= r);

        if isempty(xs_real) || isempty(rs_real), continue; end
        if length(xs_real) >= m, continue; end   % must be proper subset

        S_sub   = smat_resp(xs_real, rs_real);
        M       = length(xs_real);  N = length(rs_real);
        nullS = N - rank(full(S_sub));
        nullSp = M - rank(full(S_sub'));
        % chi     = M - N + nullS;
        chi = M - N + nullS + nullSp;

        if chi ~= 0, continue; end

        % Deduplicate
        is_dup = false;
        for k = 1:length(lbs_list)
            if isequal(sort(lbs_list{k}{1}), sort(xs_real)) && ...
               isequal(sort(lbs_list{k}{2}), sort(rs_real))
                is_dup = true;  break;
            end
        end
        if ~is_dup
            lbs_list{end+1} = {xs_real, rs_real};
        end
    end

    if isempty(lbs_list)
        fprintf('  No proper LBS found - using full network.\n');
        lbs_list = {{1:m, 1:r}};
    end

    % Sort by size (smallest first - needed for nested decomposition)
    sizes = cellfun(@(x) numel(x{1}) + numel(x{2}), lbs_list);
    [~, order] = sort(sizes);
    lbs_list = lbs_list(order);

    fprintf('  Found %d LBS candidate(s)\n', length(lbs_list));
end


function bs = expand_bs_subnetwork(smat_bin, start_reaction)
% EXPAND_BS_SUBNETWORK - Grow a subnetwork from a single seed reaction.

    [~, r] = size(smat_bin);
    rs = start_reaction;
    xs = find(any(smat_bin(:, rs), 2));

    prev_xs = [];  prev_rs = [];

    while true
        new_rs = [];
        for j = 1:r
            if ~ismember(j, rs)
                regulated_xs = find(smat_bin(:, j));
                if ~isempty(regulated_xs) && all(ismember(regulated_xs, xs))
                    new_rs = [new_rs, j]; %#ok<AGROW>
                end
            end
        end
        rs = unique([rs, new_rs]);
        xs = unique([xs; find(any(smat_bin(:, rs), 2))]);

        if isequal(prev_rs, rs) && isequal(prev_xs, xs), break; end
        prev_rs = rs;  prev_xs = xs;
    end

    bs = struct('xs', sort(xs), 'rs', sort(rs));
end


function det_struct_list = build_determinant_structures(lbs_list, model, m, r)
% BUILD_DETERMINANT_STRUCTURES - Extract Γβ shells by nested LBS peeling.
%
% Methodology: sort LBS smallest -> largest; each Γβ = LBS_k \ (all
% smaller inner LBS already processed).  The complement of all Γβ is
% added at the end.

    det_struct_list = {};

    sizes = cellfun(@(x) numel(x{1}) + numel(x{2}), lbs_list);
    [~, order] = sort(sizes, 'ascend');
    lbs_list = lbs_list(order);

    used_sp = false(1, m);
    used_rx = false(1, r);

    fprintf('  Determinant structures:\n');

    for k = 1:numel(lbs_list)
        LBS_S = lbs_list{k}{1};
        LBS_R = lbs_list{k}{2};

        Gamma_S = LBS_S(~used_sp(LBS_S));
        Gamma_R = LBS_R(~used_rx(LBS_R));

        used_sp(Gamma_S) = true;
        used_rx(Gamma_R) = true;

        if isempty(Gamma_S) && isempty(Gamma_R), continue; end

        idx = numel(det_struct_list) + 1;
        fprintf('    Γ%d (LBS %d): \n \t species {%s}, \n \t rxns %s\n', ...
            idx, k, ...
            strjoin(model.species(Gamma_S), ', '), mat2str(Gamma_R));

        det_struct_list{end+1} = struct( ...
            'species',       Gamma_S, ...
            'rxns',          Gamma_R, ...
            'label',         sprintf('Γ%d', idx), ...
            'source_lbs',    k, ...
            'is_complement', false);
    end

    % Complement: species / reactions outside all LBS
    rem_S = find(~used_sp);  rem_S = rem_S(rem_S <= m);
    rem_R = find(~used_rx);  rem_R = rem_R(rem_R <= r);

    if ~isempty(rem_S) || ~isempty(rem_R)
        fprintf('    Γ_complement: \n \t species {%s}, \n \t rxns %s\n', ...
            strjoin(model.species(rem_S), ', '), mat2str(rem_R));

        det_struct_list{end+1} = struct( ...
            'species',       rem_S, ...
            'rxns',          rem_R, ...
            'label',         'Γ_complement', ...
            'source_lbs',    [], ...
            'is_complement', true);
    end

    fprintf('  Total: %d determinant structure(s)\n', numel(det_struct_list));
end


function [bif_species, bif_rxns, bif_sources] = analyze_bifurcations( ...
    det_struct_list, A_sym, model)
% ANALYZE_BIFURCATIONS - Check each Γβ and its complement Γβ'.

    bif_species = {};  bif_rxns = {};  bif_sources = {};

    all_species = 1:size(A_sym, 2);
    all_rxns    = 1:size(A_sym, 1);

    for d = 1:numel(det_struct_list)
        Gamma_S      = det_struct_list{d}.species;
        Gamma_R      = det_struct_list{d}.rxns;
        label        = det_struct_list{d}.label;
        source_lbs   = det_struct_list{d}.source_lbs;
        is_complement = det_struct_list{d}.is_complement;

        % Analyse Γβ itself
        [found, bsp, brx, bsrc] = analyze_block( ...
            Gamma_S, Gamma_R, label, source_lbs, A_sym, model, false);
        if found
            bif_species = [bif_species, bsp];
            bif_rxns    = [bif_rxns,    brx];
            bif_sources = [bif_sources, bsrc];
        end

        % Analyse complement Γβ' (skip if this IS the complement)
        % if ~is_complement
        %     Gamma_S_comp = setdiff(all_species, Gamma_S);
        %     Gamma_R_comp = setdiff(all_rxns, Gamma_R);
        %     label_comp = [label, ' (Complement)'];
        % 
        %     [found_bif, bif_sp, bif_rx, bif_src] = analyze_block(...
        %         Gamma_S_comp, Gamma_R_comp, label_comp, source_lbs, ...
        %         A_sym, model, true);
        % 
        %     if found_bif
        %         bif_species = [bif_species, bif_sp];
        %         bif_rxns = [bif_rxns, bif_rx];
        %         bif_sources = [bif_sources, bif_src];
        %     end
        % end
    end

    if isempty(bif_species)
        fprintf('\nNo bifurcation candidates detected.\n');
    end
end


function [found, bif_sp, bif_rx, bif_src] = analyze_block( ...
    S, R, label, source_lbs, A_sym, model, is_complement_analysis)
    % ANALYZE_BLOCK - Evaluate one determinant structure for bifurcation.

    found = false;  bif_sp = {};  bif_rx = {};  bif_src = {};

    if isempty(S) || isempty(R), return; end

    [full_r, full_m] = size(A_sym);
    S_comp = setdiff(1:full_m, S);
    R_comp = setdiff(1:full_r, R);

    % Rearrange A so that Γβ occupies the top-left corner
    row_perm = [R, R_comp];
    col_perm = [S, S_comp];
    A_rearranged = A_sym(row_perm, col_perm);

    nr = length(R);
    nc = length(S);

    fprintf('\n[Analysis] %s: %d rxns × %d species\n', label, nr, nc);

    if nr == nc
        fprintf('  Square block - det(A_Γβ)\n');
        A_block = A_rearranged(1:nr, 1:nc);
    else
        fprintf('  Rectangular (%d×%d) - analysing full rearranged matrix\n', nr, nc);
        A_block = A_rearranged;
    end

    detA = simplify(det(A_block));

    if detA == 0, return; end

    factors = factor(detA);
    A_gamma_block = A_rearranged(1:nr, 1:nc);

    for k = 1:numel(factors)
        vars = symvar(factors(k));

        gamma_vars = symvar(A_gamma_block);
        factor_vars_in_gamma = intersect(vars, gamma_vars);

        if isempty(factor_vars_in_gamma) && nr ~= nc, continue; end

        if ~check_sign_change_possible(factors(k), vars), continue; end

        % == Bifurcation found ==========================================
        fprintf('\n┌=====================================┐\n');
        fprintf('│  Bifurcation Candidate Found!       │\n');
        fprintf('└=====================================┘\n');
        fprintf('  Structure : %s\n', label);
        fprintf('  Condition : %s = 0\n', char(factors(k)));

        bif_cols = [];
        for v = 1:numel(factor_vars_in_gamma)
            [~, col] = find(A_gamma_block == factor_vars_in_gamma(v) | ...
                            A_gamma_block == -factor_vars_in_gamma(v));
            bif_cols = [bif_cols, col(:)']; %#ok<AGROW>
        end
        bif_cols      = unique(bif_cols);
        bif_sp_idx    = S(bif_cols(bif_cols <= nc));

        fprintf('  Species   : %s\n', strjoin(model.species(bif_sp_idx), ', '));
        fprintf('  Reactions : r%s\n', sprintf('%d ', R));

        found = true;
        bif_sp{end+1}  = bif_sp_idx;
        bif_rx{end+1}  = R;
        bif_src{end+1} = struct( ...
            'lbs_idx',         source_lbs, ...
            'det_struct',      label, ...
            'condition',       factors(k), ...
            'was_rectangular', (nr ~= nc), ...
            'is_complement',   is_complement_analysis);
    end
end


function can_change = check_sign_change_possible(factor_expr, vars)
% CHECK_SIGN_CHANGE_POSSIBLE - Can this symbolic factor equal zero?
%
% Rules (all parameters strictly positive):
%   constant         -> check value
%   single variable  -> always > 0, cannot be zero
%   all-positive sum -> always > 0, cannot be zero
%   all-negative sum -> always < 0, cannot be zero
%   mixed signs      -> can balance to zero

    if isempty(vars)
        try
            val = double(factor_expr);
            can_change = abs(val) < 1e-10;
        catch
            can_change = false;
        end
        return;
    end

    if isscalar(vars)
        fprintf('    [Filter] Single var %s (always > 0)\n', char(factor_expr));
        can_change = false;
        return;
    end

    try
        expr_str = strrep(char(expand(factor_expr)), ' ', '');

        plus_pos  = strfind(expr_str, '+');
        minus_pos = strfind(expr_str, '-');

        has_lead_minus    = ~isempty(minus_pos) && minus_pos(1) == 1;
        n_internal_minus  = length(minus_pos) - has_lead_minus;
        n_plus            = length(plus_pos);

        if n_plus > 0 && n_internal_minus > 0
            fprintf('    [Check] Possible sign change %s\n', char(factor_expr));
            can_change = true;   % mixed signs
        elseif n_plus == 0 && n_internal_minus == 0
            fprintf('    [Filter] No operators: %s\n', char(factor_expr));
            can_change = false;  % single product, always > 0
        elseif n_plus > 0 && n_internal_minus == 0
            fprintf('    [Filter] All-positive: %s\n', char(factor_expr));
            can_change = false;
        elseif has_lead_minus && n_plus == 0
            fprintf('    [Filter] All-negative: %s\n', char(factor_expr));
            can_change = false;
        else
            fprintf('    [Check] Possible sign change %s\n', char(factor_expr));
            can_change = true;   % e.g. "r3 - r5 - r7"
        end
    catch
        can_change = true;   % conservative fallback
    end
end

%% ===================================================================
%  SYMBOLIC UTILITIES
% ===================================================================
function A_sym = symbolify_amat(amat, model, r, m)
% SYMBOLIFY_AMAT - Replace each nonzero entry amat(i,j) with a named
% symbolic variable r{i}_{species_j}, preserving the sign.

    A_sym = sym(amat);
    for i = 1:r
        for j = 1:m
            if amat(i, j) ~= 0
                % sym_name = sym(sprintf('r%d', i));
                sym_name = sym(sprintf('r%d_%s', i, model.species{j}));
                A_sym(i, j) = sign(amat(i, j)) * sym_name;
            end
        end
    end
end

function out = randomize_nonzero_entries(mat)
    % RANDOMIZE_NONZERO_ENTRIES: Generic positioning
    
    mask = (mat ~= 0);
    rand_mag = 1 + 0.5 * rand(size(mat));
    rand_signed = sign(mat) .* rand_mag;
    
    out = mat;
    out(mask) = rand_signed(mask);
end

%% ===================================================================
%  DISPLAY UTILITIES
% ===================================================================
function display_lbs(lbs_list, model)
    fprintf('\n┌========================================┐\n');
    fprintf('│  Labeled Buffering Structures (LBS)    │\n');
    fprintf('└========================================┘\n\n');
    for k = 1:length(lbs_list)
        xs = lbs_list{k}{1};  rs = lbs_list{k}{2};
        fprintf('LBS %d:\n', k);
        fprintf('  Species   : {%s}\n', strjoin(model.species(xs), ', '));
        fprintf('  Reactions : %s\n\n', mat2str(rs));
    end
end


function display_bifurcation_summary(bif_species, bif_rxns, bif_sources, model, lbs_list)
    if isempty(bif_species)
        fprintf('\n╔════════════════════════════════════════╗\n');
        fprintf('║  No bifurcations detected              ║\n');
        fprintf('╚════════════════════════════════════════╝\n');
        return;
    end

    % Deduplicate by condition string
    seen_conditions = {};
    unique_bifs = struct('condition', {}, 'species', {}, 'rxns', {}, ...
                         'det_struct', {}, 'lbs_idx', {});

    for i = 1:numel(bif_species)
        cond = char(bif_sources{i}.condition);
        if ~ismember(cond, seen_conditions)
            seen_conditions{end+1} = cond;

            min_lbs = find_minimal_lbs(bif_species{i}, lbs_list);
            nb = numel(unique_bifs) + 1;
            unique_bifs(nb).condition  = cond;
            unique_bifs(nb).species    = bif_species{i};
            unique_bifs(nb).rxns       = bif_rxns{i};
            unique_bifs(nb).det_struct = bif_sources{i}.det_struct;
            unique_bifs(nb).lbs_idx    = min_lbs;
        end
    end

    fprintf('\n╔════════════════════════════════════════╗\n');
    fprintf('║  BIFURCATION ANALYSIS SUMMARY          ║\n');
    fprintf('╚════════════════════════════════════════╝\n\n');
    fprintf('%d unique bifurcation condition(s) found\n\n', numel(unique_bifs));

    for i = 1:numel(unique_bifs)
        fprintf('Bifurcation %d\n', i);
        fprintf('  Condition : %s = 0\n', unique_bifs(i).condition);
        fprintf('  Species   : %s\n',     strjoin(model.species(unique_bifs(i).species), ', '));
        fprintf('  Reactions : r%s\n',    sprintf('%d ', unique_bifs(i).rxns));
        if ~isempty(unique_bifs(i).lbs_idx)
            lbs = lbs_list{unique_bifs(i).lbs_idx};
            fprintf('  Min LBS #%d: {%s}, rxns %s\n', ...
                    unique_bifs(i).lbs_idx, ...
                    strjoin(model.species(lbs{1}), ', '), mat2str(lbs{2}));
        else
            fprintf('  No Minimal LBS containing all species found.')
        end
        fprintf('\n');
    end
    fprintf('═════════════════════════════════════════\n');
end


function min_idx = find_minimal_lbs(bif_sp, lbs_list)
    min_idx = [];  min_size = inf;
    for i = 1:numel(lbs_list)
        if all(ismember(bif_sp, lbs_list{i}{1}))
            sz = numel(lbs_list{i}{1});
            if sz < min_size
                min_size = sz;  min_idx = i;
            end
        end
    end
end


function inspect_model(model)
    fprintf('\n╔════════════════════════════════════════════════════╗\n');
    fprintf('║              MODEL INSPECTION REPORT               ║\n');
    fprintf('╚════════════════════════════════════════════════════╝\n\n');
    fprintf('Model ID : %s\n\n', model.id);

    [model, m] = modelSpecies(model);
    [~, ~, ~, ~, r] = stoichMatrix(model, m);

    fprintf('┌=====================================┐\n');
    fprintf('│  SPECIES  (%d total)                 │\n', m);
    fprintf('└=====================================┘\n');
    for i = 1:m
        fprintf('  %2d.  %s\n', i, model.species{i});
    end

    fprintf('\n┌==================================================┐\n');
    fprintf('│  REACTIONS  (%d total)                            │\n', r);
    fprintf('└==================================================┘\n');
    rxn_idx = 1;
    for i = 1:numel(model.reaction)
        fprintf('  r%-4d  %s\n', rxn_idx, model.reaction(i).id);
        rxn_idx = rxn_idx + 1;
        if model.reaction(i).reversible
            fprintf('  r%-4d  %s_rev\n', rxn_idx, model.reaction(i).id);
            rxn_idx = rxn_idx + 1;
        end
    end
    fprintf('\n');
end


%% ===================================================================
%  MODEL CONSTRUCTION
% ===================================================================
function [model, m] = modelSpecies(model)
    sp = {};
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            sp{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            sp{end+1} = model.reaction(i).product(j).species;
        end
    end
    model.species = unique(sp);
    m = numel(model.species);
end


function [N, kinetic_N, reactant_complex, product_complex, r] = stoichMatrix(model, m)
    N = [];  kinetic_N = [];
    reactant_complex = [];  product_complex = [];

    for i = 1:numel(model.reaction)
        rc = zeros(m, 1);  krc = zeros(m, 1);
        for j = 1:numel(model.reaction(i).reactant)
            idx = find(strcmp(model.reaction(i).reactant(j).species, model.species), 1);
            rc(idx)  = model.reaction(i).reactant(j).stoichiometry;
            krc(idx) = model.reaction(i).kinetic.reactant1(j);
        end

        pc = zeros(m, 1);  kpc = zeros(m, 1);
        for j = 1:numel(model.reaction(i).product)
            idx = find(strcmp(model.reaction(i).product(j).species, model.species), 1);
            pc(idx)  = model.reaction(i).product(j).stoichiometry;
            kpc(idx) = model.reaction(i).kinetic.reactant2(j);
        end

        reactant_complex(:, end+1) = rc;
        product_complex(:, end+1)  = pc;
        N(:, end+1)       = pc - rc;
        kinetic_N(:, end+1) = krc;

        if model.reaction(i).reversible
            reactant_complex(:, end+1) = pc;
            product_complex(:, end+1)  = rc;
            N(:, end+1)         = rc - pc;
            kinetic_N(:, end+1) = kpc;
        end
    end

    r = size(N, 2);
end
