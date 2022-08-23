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

% 生成随机可行解
% 从可行的一个解开始
P = eye(K);
% V_RF = ones(N,N_RF);
tt = 2*pi*rand(1,N*N_RF);
ttt = exp(1j*tt); 
V_RF = reshape(ttt,N,N_RF);
% 生成随机数，才能保证V_RF满秩，否则后面的A_j不满秩，没法取逆
% N_loop = 150;
% N_loop_ = 3;
Sum_R = zeros(1, 11);
 
for SNR_dB = 4:2:10
    SNR = 10^(SNR_dB/10);
    sigma2 = P_t/SNR;
    while 1 % 判断第二次收敛（改变功率分配）
%         V_RF_ll = V_RF; % 存储上一次V_RF，用于判断收敛
%     for Nloop_ = 1:N_loop_
        H_t = (P)^(-0.5) * H; 
        temp_R2 = 0;
        while 1 % 判断第一次收敛
%         for Nloop = 1:N_loop
%             V_RF_last = V_RF; % 存储上一次V_RF，用于判断收敛
%             V_RF = change_V_RF(beta,K,P,sigma2); % 更新一次V_RF
            for j = 1:N_RF
                % 计算A_j
                V_RF_bar = V_RF;
                V_RF_bar(:,j) = [];
                A{j} = H_t * V_RF_bar*(V_RF_bar')* H_t';
                % 计算B_j,D_j,用于计算zeta和eta
                A_inv = (A{j})^(-1);
%                 A_inv = pinv(A{j});
                B{j} = H_t'* A_inv^2 * H_t;
                D{j} = H_t'* A_inv * H_t;
                for i = 1:N
                    % 计算zeta和eta
            %         temp = 0;
            %         for m = 1:1:N
            %             if( m~=i )
            %                 for n = 1:1:N
            %                     if(n~=i)
            %                         temp = temp + conj(V_RF(m,j))*B{j}(m,n)*V_RF(n,j);
            %                     end
            %                 end
            %             end
            %         end
                    % 计算zeta_B
                    V_RF_i = V_RF;
                    V_RF_i(i,:) = []; % 去掉第i行
                    B_temp = B{j};
                    B_temp(i,:) = [];
                    B_temp(:,i) = []; % 去掉第i行，第i列
                    VV_B = V_RF_i(:,j)' * B_temp * V_RF_i(:,j);
                    zeta_B(i,j)= B{j}(i,i) + 2 * real(VV_B);
            %         temp = 0;
            %         for m = 1:1:N
            %             if( m~=i )
            %                 for n = 1:1:N
            %                     if(n~=i)
            %                         temp = temp + conj(V_RF(m,j))*D{j}(m,n)*V_RF(n,j);
            %                     end
            %                 end
            %             end
            %         end
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
%             flag = 1;
%             for jj = 1:1:N_RF
%                 for ii = 1:1:N
%                     distance = abs( V_RF(ii,jj) - V_RF_last(ii,jj) )/abs(V_RF_last(ii,jj));
%                     if(distance >= 0.05) 
%                         flag = 0; 
%                         break;
%                     end
%                 end
%                 if(flag == 0) break; end
%             end
%             if(flag == 1) 
%                 disp("successfully converged"); 
%                 V_RF;
%                 break; 
%             end
            
        end

    %     max_ = 0;
    %     for jj = 1:N_RF
    %         for ii = 1:N
    %             ttt_r = abs(real(V_RF(ii,jj))/real(V_RF_last(ii,jj))-1);
    %             ttt_i = abs(imag(V_RF(ii,jj))/imag(V_RF_last(ii,jj))-1);
    %             max_ = max([ttt_r,ttt_i,max_]);
    %         end
    %     end
        % 生成功率分配矩阵
        V_D_t = (V_RF') * (H')/ ( H * V_RF * (V_RF') * (H'));
        Q_t = (V_D_t') * (V_RF') * V_RF * V_D_t;

%         % 迭代求出lamda
        [lambda,temp_p] = get_lambda(Q_t,P_t,sigma2,K,beta);
%         lambda = 1/sigma2;
%         while 1
%             initPower = 0;
%             posi = 0;
%             for k = 1:1:K
%                 temp_p(k) = (beta(k)/lambda) - Q_t(k,k)*sigma2;
%                 if(temp_p(k) > 0)
%                     initPower = initPower + temp_p(k);
%                     posi = posi + 1;
%                 end
%             end
% %             initPower / P_t -1;
% %             initPower;
%             if( abs(initPower / P_t -1) <= 0.001 )
%                 disp("found P");
%                 break;
%             end
% %             if(posi > 0) lamda = lamda + 0.01*(initPower - P_t)/posi;
% %             else lamda = lamda/5;
% %             end
%             lambda = lambda + 0.001*(initPower - P_t);
%             if(lambda <= 0) 
%                 pause;
%             end
%         end
        
        % 求出P
        P = zeros(K,K);
        for k = 1:1:K
%             t = (beta(k)/lambda) - Q_t(k,k)*sigma2;
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
%         flag_2 = 1;
%         for jj = 1:1:N_RF
%             for ii = 1:1:N
%                 distance = abs( V_RF(ii,jj) - V_RF_ll(ii,jj) )/abs(V_RF_ll(ii,jj));
%                 if(distance >= 0.05) % 如果某一个值和上一次循环相差超过5%，认为不收敛，跳出循环
%                     flag_2 = 0; 
%                     break;
%                 end
%             end
%             if(flag_2 == 0) 
%                 break; 
%             end
%         end
%         if(flag_2 == 1) 
%             disp("successfully converged again"); 
%             V_RF;
%             break; 
%         end
%         Dist = (V_RF - V_RF_ll)./V_RF_ll;
%         if( abs(max(Dist)) <= 0.05)
%             break;
%         end
    end
    
    % 计算V_D
    V_D = (V_RF') * (H') / ( H *V_RF * (V_RF') *(H'));

    % 计算 R_k 求和
    R = zeros(1,K);
    for k = 1:1:K
        R(k) = beta(k) * log2(1+(P(k,k)/sigma2));
    end
    
    Sum_R(SNR_dB/2 + 6) = sum(R) % Rate
end

SNR_ = -10:2:10;
plot( SNR_ , Sum_R, '-o');