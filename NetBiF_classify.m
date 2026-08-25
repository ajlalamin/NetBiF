%% NetBiF_classify.m
% ─────────────────────────────────────────────────────────────────────────────
% Bifurcation-type classification via the flux-augmented criteria of:
%
%   Okada, A., Mochizuki, A., Furuta, M., & Tsai, J.-C.
%   Flux-augmented bifurcation analysis in chemical reaction network systems.
%   Phys. Rev. E, 103, 062212 (2021).  https://doi.org/10.1103/PhysRevE.103.062212
%
%   Explicit formulations for the case of power law kinetics are derived
%   in the paper by Alamin, Hernandez (2026).
% ─────────────────────────────────────────────────────────────────────────────
% Last Modified: 26 August 2026

function NetBiF_classify(bif_species, bif_rxns, bif_sources, ...
        lbs_list, model, amat, A_sym, smat, kmat)

    type_bif = cell(1, numel(bif_species));

    if isempty(bif_species)
        fprintf('\n  No bifurcation candidates to classify.\n');
        return;
    end

    % ── Header ───────────────────────────────────────────────────────────────
    fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
    fprintf('║  BIFURCATION TYPE CLASSIFIER                               ║\n');
    fprintf('╚════════════════════════════════════════════════════════════╝\n\n');

    fprintf('Detected bifurcation candidates:\n\n');
    for b = 1:numel(bif_species)
        src = bif_sources{b};
        lbl = 'full network / complement';
        if ~isempty(src.lbs_idx), lbl = sprintf('LBS %d', src.lbs_idx); end
        fprintf('  [%d]  %-25s  (%s)\n',  b, src.det_struct, lbl);
        fprintf('       Condition : %s = 0\n',  char(src.condition));
        fprintf('       Species   : {%s}\n',    strjoin(model.species(bif_species{b}), ', '));
        fprintf('       Reactions : r%s\n\n',   sprintf('%d ', bif_rxns{b}));
    end

    % ── User input (1): candidate ─────────────────────────────────────────────
    target = input('Select candidate to classify (0 to skip): ');
    if isempty(target) || target < 1 || target > numel(bif_species)
        fprintf('  Skipping.\n');  return;
    end

    % ── User input (2): reactions containing k_n ─────────────────────────────
    fprintf(['\nSpecify ALL reactions whose rate contains the bifurcation\n', ...
             'parameter k_n (global reaction indices).\n', ...
             'Example: α in both r_2 and r_4 → enter [2 4].\n', ...
             'Reactions in selected candidate: r%s\n'], ...
             sprintf('%d ', bif_rxns{target}));
    param_rxns = input('  Reactions containing k_n: ');
    param_rxns = param_rxns(:)';

    if isempty(param_rxns)
        fprintf('  ERROR: at least one reaction index required.\n');
        type_bif{target} = 'UNCLASSIFIED (no parameter reactions)';  return;
    end

    % ══════════════════════════════════════════════════════════════════════════
    %  STEP 0 — LBS restriction and manifold enforcement
    % ══════════════════════════════════════════════════════════════════════════
    fprintf('──────────────────────────────────────────────────────────\n');  
    fprintf(' STEP 0 : LBS restriction and manifold enforcement\n');  
    fprintf('──────────────────────────────────────────────────────────\n');

    src     = bif_sources{target};
    lbs_idx = src.lbs_idx;

    % condition is stored as a symbolic object directly — no str2sym needed
    % amat_bif = enforce_bifurcation_manifold(amat, A_sym, src.condition);
    % 
    if isempty(lbs_idx)
        fprintf('  No minimal LBS — using full network.\n');
        S_LBS    = (1:size(smat,1))';
        R_LBS    = (1:size(smat,2))';
        % A_LBS    = amat;
        kmat_LBS = kmat;
        smat_L = smat;
    else
        S_LBS = lbs_list{lbs_idx}{1}(:);
        R_LBS = lbs_list{lbs_idx}{2}(:);
        % m_l   = length(S_LBS);   r_l = length(R_LBS);
        % 
        smat_L = smat(S_LBS, R_LBS);
        % cmat_L = null(smat_L,  'rational')';
        % dmat_L = null(smat_L', 'rational')';
        % 
        % A_LBS = zeros(r_l + size(dmat_L,1),  m_l + size(cmat_L,1));
        % A_LBS(1:r_l,     1:m_l)     =  amat(R_LBS, S_LBS);
        % A_LBS(1:r_l,     m_l+1:end) =  cmat_L';
        % A_LBS(r_l+1:end, 1:m_l)     = -dmat_L;

        kmat_LBS = kmat(S_LBS, R_LBS);   % m_l × r_l,  kmat_LBS(j,i) = α_ij

        fprintf('  LBS %d : %d species  {%s}\n', ...
                lbs_idx, length(S_LBS), strjoin(model.species(S_LBS), ', '));
        fprintf('          %d reactions  r%s\n', length(R_LBS), sprintf('%d ', R_LBS));
    end

    m_lbs = length(S_LBS);
    r_lbs = length(R_LBS);

    % ── Symbolic steady state ─────────────────────────────────────────────────
    %   x̄ is kept fully symbolic — no numerical value assumed.
    x_bar_lbs_sym = sym(zeros(m_lbs, 1));
    for j = 1:m_lbs
        raw  = model.species{S_LBS(j)};
        safe = regexprep(raw, '[^a-zA-Z0-9]', '_');   % MATLAB sym-safe name
        x_bar_lbs_sym(j) = sym(sprintf('s_%s', safe), 'positive');
    end
 
    fprintf('\n  Symbolic steady state x̄ (LBS species):\n');
    for j = 1:m_lbs
        fprintf('    x̄_%s  =  %s\n', model.species{S_LBS(j)}, ...
                char(x_bar_lbs_sym(j)));
    end

    % Translate global reaction indices to local (within R_LBS)
    param_rxns_local = [];
    for k = 1:numel(param_rxns)
        loc = find(R_LBS == param_rxns(k));
        if ~isempty(loc)
            param_rxns_local(end+1) = loc; %#ok<AGROW>
        else
            fprintf('  WARNING: r%d not in LBS — ignored.\n', param_rxns(k));
        end
    end

    if isempty(param_rxns_local)
        fprintf('  ERROR: no specified reactions are inside the LBS.\n');
        type_bif{target} = 'UNCLASSIFIED (parameter reactions outside LBS)';  return;
    end

    % ══════════════════════════════════════════════════════════════════════════
    %  STEP 1 — SN1/T1/P1 : simple zero eigenvalue of A_LBS
    % ══════════════════════════════════════════════════════════════════════════
    fprintf('──────────────────────────────────────────────────────────\n');  
    fprintf(' STEP 1: Simple zero eigenvalue of A_LBS\n');  
    fprintf('──────────────────────────────────────────────────────────\n');

    % 1. Restrict A_sym to LBS
    A_sym_L = A_sym(R_LBS, S_LBS);
    cmat_L  = null(smat_L,  'rational')';
    dmat_L  = null(smat_L', 'rational')';

    % Normalise each row (each basis vector) to remove large integer factors
    for i = 1:size(cmat_L, 1)
        nz = find(cmat_L(i,:) ~= 0, 1);
        if ~isempty(nz), cmat_L(i,:) = cmat_L(i,:) / cmat_L(i,nz); end
    end
    for i = 1:size(dmat_L, 1)
        nz = find(dmat_L(i,:) ~= 0, 1);
        if ~isempty(nz), dmat_L(i,:) = dmat_L(i,:) / dmat_L(i,nz); end
    end
    
    A_sym_LBS = [A_sym_L,          sym(cmat_L'); ...
                 sym(-dmat_L),      sym(zeros(size(dmat_L,1), size(cmat_L,1)))];

    % 2. Enforce bifurcation condition symbolically
    %    (substitute target_var = f(other_vars) directly into A_sym_LBS)
    vars         = symvar(src.condition);
    target_var   = vars(1);
    sol = solve(src.condition == 0, target_var);
    A_sym_bif    = subs(A_sym_LBS, target_var, sol(1));

    A_sym_bif = simplify(A_sym_bif);

    remaining_vars = symvar(A_sym_bif);
    assume(remaining_vars, 'real');
    assume(remaining_vars, 'positive');
    
    % 3. Compute symbolic eigenvectors
    W_sym = null(A_sym_bif');   % left  eigenvector (symbolic)
    U_sym = null(A_sym_bif);    % right eigenvector (symbolic)

    % W_sym = simplify(W_sym / W_sym(find(W_sym ~= 0, 1)));
    % Filter out zero entries before selecting normalization pivot
    nonzero_mask_W = arrayfun(@(x) ~isequal(x, sym(0)), W_sym);  % element-wise
    valid_entries_W = W_sym(nonzero_mask_W);
    [~, idx] = min(cellfun(@(x) length(char(x)), num2cell(valid_entries_W)));
    W_sym = W_sym / valid_entries_W(idx);

    % U_sym = simplify(U_sym / U_sym(find(U_sym ~= 0, 1)));
    nonzero_mask_U = arrayfun(@(x) ~isequal(x, sym(0)), U_sym);
    valid_entries_U = U_sym(nonzero_mask_U);
    [~, idx] = min(cellfun(@(x) length(char(x)), num2cell(valid_entries_U)));
    U_sym = U_sym / valid_entries_U(idx);

    fprintf('  W (symbolic): %s\n', char(W_sym'));
    fprintf('  U (symbolic): %s\n', char(U_sym'));
    fprintf('──────────────────────────────────────────────────────────\n');
    
    % 4. Evaluate derivative terms symbolically
    u_sym       = U_sym(1:m_lbs);          % species block
    u_tilde = simplify(u_sym ./ x_bar_lbs_sym);   % ũ, length m_lbs
    % x_bar_sym   = sym(x_bar(1:m_lbs));

    % ensure that system is at bifurcation point
    assume(src.condition == 0);
    
    if simplify(src.condition) == 0
        at_bif = "TRUE";
    else
        at_bif = "FALSE";
    end
    fprintf('Assume dynamics at bifurcation point: %s \n', at_bif);
    
    %% PARAMETER CALCULATIONS

    Theta_sym = compute_Theta_condition(W_sym, x_bar_lbs_sym, ...
                                        kmat_LBS, param_rxns_local);

    fprintf('──────────────────────────────────────────────────────────\n');
    fprintf(' (A0) : Net parameter sensitivity  Θ\n');
    fprintf('  Θ (symbolic) = %s\n', char(Theta_sym));
    [N_theta, ~] = numden(Theta_sym);
    Theta_factors = factor(N_theta);
    % fprintf('  Θ (symbolic) = %s\n', char(N_theta));
    % Theta_factors = factor(Theta_sym);
    theta_can_vanish = any(arrayfun( ...
        @(f) check_sign_change_possible(f, symvar(f)), Theta_factors));
    if theta_can_vanish
        fprintf("     [Check] Θ CAN be zero\n");
    else
        fprintf("     [Check] Θ is structurally NONZERO\n")
    end
    fprintf('──────────────────────────────────────────────────────────\n');

    A1_sym = compute_A1_condition(W_sym, u_tilde, x_bar_lbs_sym, ...
                                  kmat_LBS, param_rxns_local);
    A1_sym = simplify(A1_sym);
 
    fprintf('──────────────────────────────────────────────────────────\n');
    fprintf(' (A1) : Mixed partial  W_n [ ∇_x(∂r_n/∂k_n) · u ]\n');
    fprintf('  A1 (symbolic) = %s\n', char(A1_sym));
    [N_A1, ~] = numden(A1_sym);
    A1_factors = factor(N_A1);
    a1_can_vanish = any(arrayfun( ...
        @(f) check_sign_change_possible(f, symvar(f)), A1_factors));
    if a1_can_vanish
        fprintf('     [Check] A1 CAN be zero\n');
    else
        fprintf('     [Check] A1 is structurally NONZERO\n');
    end
    fprintf('──────────────────────────────────────────────────────────\n');

    H_sym = compute_hessian_contraction(W_sym, u_tilde, x_bar_lbs_sym, A_sym_bif, kmat_LBS, r_lbs, m_lbs);
    H_sym = simplify(H_sym);
            
    fprintf('──────────────────────────────────────────────────────────\n');
    fprintf(' (A2) : Hessian Contraction at (u,u)\n');
    fprintf('  W^T [D²_x r(u,u)]  =  %s\n', char(H_sym));
    [N_H, ~] = numden(H_sym);
    H_factors = factor(N_H);
    % H_factors = factor(H_sym);
    h_can_vanish = any(arrayfun( ...
        @(f) check_sign_change_possible(f, symvar(f)), H_factors));
    if h_can_vanish
        fprintf("     [Check] Hessian CAN be zero\n"); 
        % S = solve(N_H,symvar(H_factors));
        % disp(S);
    else
        fprintf("     [Check] Hessian is structurally NONZER\n")
    end
    fprintf('──────────────────────────────────────────────────────────\n');

    T3_sym = compute_third_derivative_contraction(W_sym, u_tilde, x_bar_lbs_sym, A_sym_bif, kmat_LBS, r_lbs, m_lbs);
    T3_sym = simplify(T3_sym);

    fprintf('──────────────────────────────────────────────────────────\n');
    fprintf(' (A3) : Third Derivative Contraction at (u,u,u)\n');
    fprintf('  W^T [D³_x r(u,u,u)]  =  %s\n', T3_sym);
    [N_T3, ~] = numden(T3_sym);
    T3_factors = factor(N_T3);
    % T3_factors = factor(T3_sym);
        t3_can_vanish = any(arrayfun( ...
            @(f) check_sign_change_possible(f, symvar(f)), T3_factors));
    if t3_can_vanish
        fprintf("     [Check] Third Derivative CAN be zero\n"); 
    else
        fprintf("     [Check] Third Derivative is structurally NONZERO\n"); 
    end
    fprintf('──────────────────────────────────────────────────────────\n\n');

    fprintf('──────────────────────────────────────────────────────────\n');

    %% CLASSIFICATION PROPER

    if ~theta_can_vanish
        % Θ ≠ 0 structurally → SADDLE-NODE branch
        fprintf(' Θ ≠ 0 structurally → SADDLE-NODE branch\n');
        fprintf('──────────────────────────────────────────────────────────\n');
    
        if ~h_can_vanish
            bif_type = 'SADDLE-NODE';
        else
            bif_type = 'SADDLE-NODE (possible degenerate/cusp)';
        end
    
    else
        % Θ = 0 reachable → TRANSCRITICAL / PITCHFORK branch
        fprintf(' Θ = 0 reachable → TRANSCRITICAL / PITCHFORK branch\n');
        fprintf('──────────────────────────────────────────────────────────\n');

        if a1_can_vanish
            fprintf([' A1 = 0 reachable  →  condition (A1) fails\n', ...
                     ' Standard TC/PF classification does not apply.\n']);
            fprintf('──────────────────────────────────────────────────────────\n');
            
            bif_type = 'DEGENERATE (A1 = 0)';
 
        else
    
            if ~h_can_vanish
                % A2 ≠ 0 → TRANSCRITICAL
                bif_type = 'TRANSCRITICAL';
            else
                % A2 = 0 → check A3
                fprintf(' A2 = 0 reachable → checking A3 for PITCHFORK\n');
                fprintf('──────────────────────────────────────────────────────────\n');
        
                if ~t3_can_vanish
                    bif_type = 'PITCHFORK';
                else
                    fprintf(' A3 = 0 reachable → CANNOT determine type\n');
                    bif_type = 'CANNOT DETERMINE';
                end
            end
        end
    end

    %% ── Result ───────────────────────────────────────────────────────────────
    fprintf('\n╔════════════════════════════════════════════════════════════╗\n');
    fprintf('║  CLASSIFICATION: %-42s║\n', [bif_type, '  ']);
    fprintf('╚════════════════════════════════════════════════════════════╝\n');
    fprintf('  Parameter reactions : r%s\n', sprintf('%d ', param_rxns));
    fprintf('  Steady state x̄     : symbolic  (%s)\n', ...
            strjoin(arrayfun(@char, x_bar_lbs_sym, 'UniformOutput', false), ', '));
    if ~isempty(lbs_idx)
        fprintf('  Minimal LBS         : LBS %d  {%s}\n', ...
                lbs_idx, strjoin(model.species(S_LBS), ', '));
    end
    fprintf('\n');
    type_bif{target} = bif_type;

    % do_hopf = input('Run supplementary Hopf check on A_LBS? (1=yes / 0=no): ');
    % if ~isempty(do_hopf) && do_hopf == 1
    %     hopf_check(A_LBS, TOL);
    % end
end


%% =============================================================================
%%  COMPUTATION FUNCTIONS
%% =============================================================================


function Theta_sym = compute_Theta_condition(W_sym, x_bar_lbs_sym, ...
                                            kmat_LBS, param_rxns_local)
% Θ — net parameter sensitivity
    Theta_sym = sym(0);

    for i = param_rxns_local(:)'
        alpha_i   = sym(kmat_LBS(:, i));               % m_lbs × 1
        dr_dk_i   = prod(x_bar_lbs_sym .^ alpha_i);   % ∏ x̄_j^{α_ij}
        Theta_sym = Theta_sym + W_sym(i) * dr_dk_i;
    end
    
    Theta_sym = simplify(Theta_sym);
end
function A1_sym = compute_A1_condition(...
        W_sym, u_tilde, x_bar_sym, kmat_LBS, param_rxns_local)
% COMPUTE_A1_CONDITION

    A1_sym = sym(0);

    for n = param_rxns_local(:).'
        alpha_n = sym(kmat_LBS(:, n));           % m_lbs × 1

        % φ_n(x̄) = ∏_j x̄_j^{α_nj}
        phi_n = prod(x_bar_sym .^ alpha_n);      % scalar

        A1_sym = A1_sym + W_sym(n) * phi_n * (alpha_n.' * u_tilde);
    end

    A1_sym = simplify(A1_sym);
end

function H_sym = compute_hessian_contraction( ...
        W_sym, u_tilde, x_bar_sym, A_sym_bif, kmat_LBS, r_lbs, m_lbs)
% COMPUTE_HESSIAN_CONTRACTION
%
%   Computes  W^T [D²_x r(u,u)]  under power-law kinetics.
%
%   Under  r_i(x) = k_i · ∏_j x_j^{α_ij}:
%
%     ∂r_i/∂x_j |_{x̄}  =  r_i(x̄) · α_ij / x̄_j
%     ⟹  r_i(x̄)  =  (∂r_i/∂x_j)|_{x̄} · x̄_j / α_ij
%
%     D²r_i(u,u)  =  r_i(x̄) · [ (ũᵀ α_i)² − Σ_j α_ij ũ_j² ]
 
    alpha = sym(kmat_LBS);               % m_lbs × r_lbs
    drdx  = A_sym_bif(1:r_lbs, 1:m_lbs);
 
    % Recover r_i(x̄) from any nonzero-kinetic-order column:
    %   r_i(x̄) = drdx(i,j) · x̄_j / α_ij
    r_sym = sym(zeros(r_lbs, 1));
    for i = 1:r_lbs
        j = find(kmat_LBS(:, i) ~= 0, 1);
        if ~isempty(j)
            r_sym(i) = drdx(i, j) * x_bar_sym(j) / alpha(j, i);
        end
    end
 
    % α_i^T ũ  and  Σ_j α_ij ũ_j²
    Au    = alpha.' * u_tilde;          % r_lbs × 1
    Bu    = alpha.' * (u_tilde .^ 2);   % r_lbs × 1
 
    H_sym = W_sym(1:r_lbs).' * (r_sym .* (Au .^ 2 - Bu));
end
 
 
function T3_sym = compute_third_derivative_contraction( ...
        W_sym, u_tilde, x_bar_sym, A_sym_bif, kmat_LBS, r_lbs, m_lbs)
% COMPUTE_THIRD_DERIVATIVE_CONTRACTION
%
%   Computes  W^T [D³_x r(u,u,u)]  under power-law kinetics.
%
%     D³r_i(u,u,u) = r_i(x̄) · [ (ũᵀα_i)³
%                                  − 3(ũᵀα_i)(Σ_j α_ij ũ_j²)
%                                  + 2(Σ_j α_ij ũ_j³) ]
 
    alpha = sym(kmat_LBS);
    drdx  = A_sym_bif(1:r_lbs, 1:m_lbs);
 
    r_sym = sym(zeros(r_lbs, 1));
    for i = 1:r_lbs
        j = find(kmat_LBS(:, i) > 0, 1);
        if ~isempty(j)
            r_sym(i) = drdx(i, j) * x_bar_sym(j) / alpha(j, i);
        end
    end
 
    Au     = alpha.' * u_tilde;
    Bu     = alpha.' * (u_tilde .^ 2);
    Cu     = alpha.' * (u_tilde .^ 3);
    T3_sym = W_sym(1:r_lbs).' * ...
                 (r_sym .* (Au .^ 3 - 3 .* Au .* Bu + 2 .* Cu));
end

%% =============================================================================
%%  HOPF CHECK  (supplementary — not from Okada 2021)
%% =============================================================================

function can_change = check_sign_change_possible(factor_expr, vars)
% CHECK_SIGN_CHANGE_POSSIBLE — Can this symbolic factor equal zero?
%
% Rules (all parameters strictly positive):
%   constant         → check value
%   single variable  → always > 0, cannot be zero
%   all-positive sum → always > 0, cannot be zero
%   all-negative sum → always < 0, cannot be zero
%   mixed signs      → can balance to zero

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
