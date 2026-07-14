function icc_val = subfun_compute_icc3_1(X)
    % X: a (n_subjects x 2) matrix, columns are sessions
    if size(X,2) ~= 2
        error('Input must have exactly two columns (2 sessions).');
    end

    n = size(X, 1);  % number of subjects

    % Means
    mean_rows = mean(X, 2);    % per subject
    mean_cols = mean(X, 1);    % per session
    grand_mean = mean(X(:));   % grand mean

    % Sum of Squares
    SS_total = sum((X(:) - grand_mean).^2);
    SS_subject = 2 * sum((mean_rows - grand_mean).^2);
    SS_session = n * sum((mean_cols - grand_mean).^2);
    SS_error = SS_total - SS_subject - SS_session;

    % Degrees of Freedom
    df_subject = n - 1;
    df_session = 1; % since only 2 sessions
    df_error = (n - 1) * 1;

    % Mean Squares
    MS_subject = SS_subject / df_subject;
    MS_error = SS_error / df_error;

    % ICC(3,1) formula (consistency, fixed effects model)
    icc_val = (MS_subject - MS_error) / (MS_subject + MS_error);
end