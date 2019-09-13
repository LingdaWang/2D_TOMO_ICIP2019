clear
close all
clc

%% Load dataset
load clean.mat
img_orig = projs(:,:,101);
img = imresize(img_orig,65/129);% downsample the image
size_img = 65;% size of the image.
R = 32;
c = 0.3;% bandlimit
 
%% some parameters of img
n_r = ceil(4*c*R); % num of points in the projection lines
num_pool = 4;
[ basis, sample_points, ~ ] = Precomp_FB( n_r, R, c ); 
[ coeff, ~ ] = Coeff_FB(img, R, 0, basis, sample_points, num_pool);% get the ground truth a_kq of the img
k_max = size(coeff,1) - 1; % max frequency
index_coeff = zeros(1+k_max,1);% structure of a_kq
for i=1:1:k_max+1   
    index_coeff(i,1)=size(coeff{i},1);  
end
index_coeff_full = [fliplr(index_coeff'),(index_coeff(2:end))'];% full structure of a_kq
coeff_structure = cellfun(@(x,y)x./y,coeff,coeff,'UniformOutput',false);% structure of the a_kq

%% Generate the non_uniform pmf
delta_alpha = 2*pi/(2*k_max+1);% change to 100
alpha = 20/360.*2*pi;
num_of_delta_alpha = ceil(alpha./delta_alpha);
num_of_sample = 10000;
% Generate the pmfimg
[ pmf_theta, num_of_sample,rand_num] = Genarate_Nonuniform_Distribution( k_max, num_of_sample );
FFT_pmf = fft(pmf_theta); %% FFT of the pmf

%% Get the projection lines at differnet angle
[ F_img_eq, freqs2] = cryo_pft( img,2*R,2*k_max+1 ); 
img_real_line = real(ifft([F_img_eq;flipud(conj(F_img_eq(2:end,:)))]));
img_real_line = flipud([img_real_line(end-R+1:end,:);img_real_line(1:R+1,:)]);

%% Clean data Genetation
clean_data = clean_data_gen(img_real_line,rand_num, num_of_delta_alpha,k_max);

%% compute the weight

W_tmp = (sample_points.r.*sample_points.w).^0.5;
W_tmp = [flipud(W_tmp);W_tmp];
W = W_tmp*W_tmp';
W = repmat(W,2.*num_of_delta_alpha+1,2.*num_of_delta_alpha+1);

%% Compute the theoretical first and second moment
[a, ~] = a_Construction(coeff); % vectorize the coeff
[ Psi, ~ ] = Psi_Construction( basis, num_of_delta_alpha, delta_alpha, k_max );% Psi matrix
Psi_W = Psi.*repmat(W_tmp,2.*num_of_delta_alpha+1,1); % weighted Psi matrix
[ index_C_P ,C_P, ~ ] = C_P_Construction( FFT_pmf, coeff_structure);% matrix of probability
u_i = floor(size(C_P,1)/2);
Cov_t = Psi*((a*(a')).*C_P)*(Psi');% theoretical second moment
[Mu_t] = Mu_t_construction(Psi,a,C_P);% theoretical first moment
Cov_t_W = Psi_W*((a*(a')).*C_P)*(Psi_W');% weighted theoretical second moment
[Mu_t_W] = Mu_t_construction(Psi_W,a,C_P);% weighted theoretical first moment

%% Compute the empririal first and second moment

noise_var = 0; % noise added to each projection line.
noise = normrnd(0,sqrt(noise_var),2*R+1,size(clean_data,2));
check_var = var(reshape(noise,1,[]));
[ freqs ] = pft_freqs(sample_points.r, 2*k_max+1);
freqs = -2.*pi.*freqs(1:n_r,2).';
freqs = [-fliplr(freqs),freqs];
var_clean_data = var(reshape(clean_data,1,[]));
SNR = var_clean_data/noise_var;
noisy_data = clean_data+noise;
% clear clean_data noise
%% by regular way
Sigma = Sigma_constr(noise_var,R,n_r,freqs);
[Cov_E, Mu_E] = Cov_E_and_Mu_E(noisy_data,freqs,num_of_delta_alpha,n_r,Sigma);
Cov_E_W = Cov_E.*W; % weighted second moment
Mu_E_W = Mu_E.*repmat(W_tmp,2.*num_of_delta_alpha+1,1);% weighted first moment

%% sanity check of the empirical and theoretical
error1 = norm(Cov_t-Cov_E,'fro')./norm(Cov_t,'fro');
error2 = norm(Mu_t-Mu_E,'fro')./norm(Mu_t,'fro');
error3 = norm(Cov_t_W-Cov_E_W,'fro')./norm(Cov_t_W,'fro');
error4 = norm(Mu_t_W-Mu_E_W,'fro')./norm(Mu_t_W,'fro');

%% ADMM optimization process
%hyper parameter
% 1 1000 1 for 16*16 
% 1 30 1 for 32*32 
% 1 5 1 for 64*64 
% 1 0.5 1 for 128*128 beta=0
lambda_1 = 1;
lambda_2 = 5;
rho = 1;

% parameter initialization 

u_init = 1/rho.*zeros(size(a,1),1);

x_cell = cell(2*k_max+1,1);
for i=1:k_max+1
    if i==1
        x_cell{i+k_max,1} = 1.*abs(randn(index_coeff(i),1));
    else
        x_cell{i+k_max,1} = 0.1.*randn(index_coeff(i),1)+0.1*1i.*randn(index_coeff(i),1);
        x_cell{k_max+2-i,1} = conj(x_cell{i+k_max,1});
    end
end
x_init = cell2mat(x_cell);
y_init = x_init;
[ tmp_p_init, ~,~] = Genarate_Nonuniform_Distribution( k_max, num_of_sample );
p_init = fft(tmp_p_init);
p_init = p_init(2:end);
[ ~,C_P_init, ~ ] = C_P_Construction( [1;p_init], coeff_structure);

%precompute tha vectorize
Psi_Psi = Psi_W'*Psi_W;
b_x_1 = Mu_E_W;
b_x_2 = reshape(Cov_E_W,[],1);
Psi_b_Psi = Psi_W'*Cov_E_W*Psi_W;
Psi_b_x_1 = Psi_W'*b_x_1;
b_y_1 = b_x_1;
b_y_2 = b_x_2;
A_x_3 = diag(ones(size(a,1),1));

% iterations
[cost_track,cost1,cost2,cost3,cost4] = cost_func(lambda_1,lambda_2,rho,u_init,x_init,y_init,C_P_init,Mu_E_W,Cov_E_W,Psi_W,u_i);
for iter = 0:200
    if mod(iter,10)==0
        [cost_func_value,cost1,cost2,cost3,cost4] = cost_func(lambda_1,lambda_2,rho,u_init,x_init,y_init,C_P_init,Mu_E_W,Cov_E_W,Psi_W,u_i);
        cost_track = [cost_track,cost_func_value];
        fprintf('The %dth iteration\n',iter)
        fprintf('cost function = %d, %d, %d, %d, %d\n',cost_func_value,cost1,cost2,cost3,cost4);
    end
    if cost_func_value<1e-8
        break;
    end
    %% update x
    [A_x_2,b_x_2] = A_x_2_costruction(y_init,C_P_init,Psi_Psi,Psi_b_Psi);
    b_x_3 = y_init-u_init;
    A_x_1_A_x_1 = Psi_Psi.*(C_P_init(u_i,:).'*conj(C_P_init(u_i,:)));
    A_x_1_b_x_1 = Psi_b_x_1.*(C_P_init(u_i,:).');
    A_x = lambda_1.*A_x_1_A_x_1+lambda_2.*A_x_2+rho.*A_x_3;
    b_x = lambda_1.*A_x_1_b_x_1+lambda_2.*b_x_2+rho.*b_x_3;
    x_init = A_x\b_x;
    %% update y
    [A_y_2,b_y_2] = A_y_2_costruction(x_init,C_P_init,Psi_Psi,Psi_b_Psi);    
    b_y_3 = x_init+u_init;
    A_y = lambda_1.*A_x_1_A_x_1+lambda_2.*A_y_2+rho.*A_x_3;
    b_y = lambda_1.*A_x_1_b_x_1+lambda_2.*b_y_2+rho.*b_y_3;
    y_init = A_y\b_y;
    %% update p
    [A_p_1,tmp_p_1] = A_p_1_costruction(Psi_W,x_init,index_coeff_full,k_max);
    b_p_1 = b_x_1 - tmp_p_1;
    [A_p_2,tmp_p_2] = A_p_1_costruction(Psi_W,y_init,index_coeff_full,k_max);
    b_p_2 = b_x_1 - tmp_p_2;
    [A_p_3,b_p_3] = A_p_3_construction(Psi_W,x_init,y_init,index_coeff_full,k_max,Cov_E_W);
    A_p = lambda_1.*A_p_1'*A_p_1+lambda_1.*A_p_2'*A_p_2+lambda_2.*A_p_3;
    b_p = lambda_1.*A_p_1'*b_p_1+lambda_1.*A_p_2'*b_p_2+lambda_2.*b_p_3;
    p_init = A_p\b_p;
    [ ~,C_P_init, ~ ] = C_P_Construction( [1;p_init], coeff_structure);
    %% update u
    u_init = u_init+(x_init-y_init);
end
zeropad = 360000;
[estimated_coeff,estimated_pmf,hat_a] = coeff_pmf_recon(x_init,p_init,index_coeff_full,k_max,coeff,zeropad,pmf_theta);
estimated_pmf = ProjectOntoSimplex(estimated_pmf, 1);
img_ori_bessel = recon_images_FB(c, R, size_img, coeff, 1, 1); % image use true a_kq
img_recon = recon_images_FB(c, R, size_img, estimated_coeff, 1, 1);
img_recon_r = recon_images_FB(c, R, size_img, hat_a, 1, 1);% recon image with rotation
error = norm(img_ori_bessel-img_recon,'fro')/norm(img_ori_bessel,'fro');
fprintf('the error is %d\n',error);
