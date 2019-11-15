function C = C_ML_Lasso( y, S_init,lambda)

y = double(y);
S = S_init;
[height, width, wavelength] = size(y);
[~,k] = size (S_init); % k is number of pure components
C = zeros(height, width,k);
%e = nan(height, width, wavelength);
mask = nan(height, width, wavelength);
mask(~isnan(y)) = 1;
% Initialize e as y-A(S\circ I)C
%e = zeros(height,width,wavelength);

% Initialize C as the ML estimate, with lasso as regularization
C_prev = zeros(k,1);
sampled_sum = 0;
%squared_error_sum = 0;
for i = 1: height
    for j = 1:width
        non_zero_locs = find(mask(i,j,:)==1);
        non_zero_count = length(non_zero_locs);
        sampled_sum = sampled_sum + non_zero_count;
        S_temp = zeros(non_zero_count,k);
        y_temp = zeros(non_zero_count,1);
        for l = 1:non_zero_count
            S_temp(l,:) = S(non_zero_locs(l),:);
            y_temp(l) = y(i,j,non_zero_locs(l));
        end
        % Note that Lasso requires at least two rows for y_temp
        if non_zero_count >= 2
            C_temp = lasso(S_temp,y_temp,'lambda',lambda);
            %C_temp = S_temp\y_temp;
        else
            % Take the C of previous pixel as a rough initialization
            C_temp = C_prev;
        end
        C_prev = C_temp;
        %e_temp = y_temp - S_temp*C_temp;
        for l = 1:k
            C(i,j,l) = C_temp(l);
        end
        %for l = 1:non_zero_count
        %    e(i,j,non_zero_locs(l)) = e_temp(l);
        %end
        %squared_error_sum = squared_error_sum + sum(e_temp.^2);
    end
end


%Display the initial estimation of C
%display_results_fungi;

end