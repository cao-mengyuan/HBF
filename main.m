clc
clear

% global N N_RF
% global H_t
% global V_RF

K = 8; % 用户数
N = 64;
N_RF = 9;

P_t = K; % 总功率
beta = ones(1,K);% beta_k = 1; %均匀分配权重
% sigma2 = K;
% SNR = P_t/sigma2;% SNR = P/sigma^2;
% SNR_dB = 10* log10(SNR);

H = channel(N,K); % 生成信道
% H = conj(H);
% 生成随机可行解
% 从可行的一个解开始
P = eye(K);
% V_RF = ones(N,N_RF);
tt = 2*pi*rand(1,N*N_RF);
ttt = exp(1j*tt); 
V_RF = reshape(ttt,N,N_RF);
% 生成随机数，才能保证V_RF满秩，否则后面的A_j不满秩，没法取逆
Sum_R = zeros(1, 11);
 
for SNR_dB = -10:2:-10
    SNR = 10^(SNR_dB/10);
    sigma2 = P_t/SNR;
    while 1 % 判断第二次收敛（改变功率分配）
        H_t = (P)^(-0.5) * H; 
        temp_R2 = 0;
        while 1 % 判断第一次收敛
            for j = 1:N_RF
                % 计算A_j
                V_RF_bar = V_RF;
                V_RF_bar(:,j) = [];
                A{j} = H_t * V_RF_bar * (V_RF_bar')* H_t';
                % 计算B_j,D_j,用于计算zeta和eta
                A_inv = (A{j})^(-1);
                B{j} = H_t'* A_inv^2 * H_t;
                D{j} = H_t'* A_inv * H_t;
                for i = 1:N
                    % 计算zeta和eta
                    % 计算zeta_B
                    V_RF_i = V_RF;
                    V_RF_i(i,:) = []; % 去掉第i行
                    B_temp = B{j};
                    B_temp(i,:) = [];
                    B_temp(:,i) = []; % 去掉第i行，第i列
                    VV_B = V_RF_i(:,j)' * B_temp * V_RF_i(:,j);
                    zeta_B(i,j)= B{j}(i,i) + 2 * real(VV_B);
                    
                    % 计算zeta_D
                    D_temp = D{j};
                    D_temp(i,:) = [];
                    D_temp(:,i) = [];
                    VV_D = V_RF_i(:,j)' * D_temp * V_RF_i(:,j);
                    zeta_D(i,j)= D{j}(i,i) + 2 * real(VV_D);

                    % 计算eta_B, eta_D
                    V_B = B{j}(i,:) * V_RF(:,j);
                    V_D = D{j}(i,:) * V_RF(:,j);
                    eta_B(i,j) = V_B - B{j}(i,i)*V_RF(i,j);
                    eta_D(i,j) = V_D - D{j}(i,i)*V_RF(i,j);

                    % 计算theta_1，theta_2
                    c(i,j) = (1 + zeta_D(i,j)) * eta_B(i,j) - zeta_B(i,j) * eta_D(i,j);
                    z(i,j) = imag( 2 * eta_B(i,j) * eta_D(i,j));
                    phi(i,j) = 0;
                    tt = asin(imag(c(i,j))/abs(c(i,j)));
                    if(real(c(i,j))>=0) phi(i,j) = tt; 
                    else phi(i,j) = pi - tt; 
                    end
                    theta_1 = -phi(i,j) + asin(z(i,j)/abs(c(i,j)));
                    theta_2 = pi - phi(i,j) - asin(z(i,j)/abs(c(i,j)));

                    % 判断最优theta
                    V_RF1 = exp(-1j * theta_1);
                    V_RF2 = exp(-1j * theta_2);
                    trace_A_inv = trace(A_inv);
                    f1 = trace_A_inv - ...
                        ( zeta_B(i,j) + 2 * real(conj(V_RF1)*eta_B(i,j)))/...
                        (1+zeta_D(i,j)+ 2 * real(conj(V_RF1)*eta_D(i,j)));
                    f2 = trace_A_inv - ...
                        ( zeta_B(i,j) + 2 * real(conj(V_RF2)*eta_B(i,j)))/...
                        (1 + zeta_D(i,j)+ 2 * real(conj(V_RF2)*eta_D(i,j)));
                    if(f1 < f2) theta_opt = theta_1; 
                    else theta_opt = theta_2;
                    end

                    V_RF(i,j) = exp(-1j*theta_opt); % 更新V_RF
                end
            end

            % R = zeros(1,K);
            for k = 1:1:K
                R(k) = beta(k) * log2(1+(P(k,k)/sigma2));
            end
            sum_r = sum(R)
            
            % 判别收敛
            if(abs(temp_R2 / sum_r -1) <= 0.0001)
                disp("successfully converged"); 
                Sum_R(SNR_dB/2 + 6) = sum_r;
                break;
            end
            temp_R2 = sum_r;
        end

        % 生成功率分配矩阵
        V_D_t = (V_RF') * (H') / ( H * V_RF * (V_RF') * (H'));
        Q_t = (V_D_t') * (V_RF') * V_RF * V_D_t;

        % 迭代求出lambda
        [lambda,temp_p] = get_lambda(Q_t,P_t,sigma2,K,beta);
        
        % 求出P
        P = zeros(K,K);
        for k = 1:1:K
            if(temp_p(k) > 0)
                P(k,k) = temp_p(k)/Q_t(k,k);
            else
                P(k,k) = 0.001;
            end
        end
        
        for k = 1:1:K
            R(k) = beta(k) * log2(1+(P(k,k)/sigma2));
        end
        sum_r = sum(R)
        
        % 检查收敛
        if(abs(Sum_R(SNR_dB/2 + 6) / sum_r -1) <= 0.001)
            disp("successfully converged again"); 
            Sum_R(SNR_dB/2 + 6) = sum_r;
            break; 
        end
        Sum_R(SNR_dB/2 + 6) = sum_r;
    end
    
    % 计算V_D
    V_D = (V_RF') * (H') / ( H *V_RF * (V_RF') *(H')) * P^(0.5);

    % 计算 R_k 求和
    R = zeros(1,K);
    for k = 1:1:K
        R(k) = beta(k) * log2(1+(P(k,k)/sigma2));
    end
    
    Sum_R(SNR_dB/2 + 6) = sum(R) % Rate
end

SNR_ = -10:2:10;
plot( SNR_ , Sum_R, '-o');