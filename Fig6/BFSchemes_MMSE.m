function [MMSE_rates,MMSE]=BFSchemes_MMSE(Channel_mul,var_noise,power,user_num,L,N,M)
 
%% %ZF for Perfect channels
%%预编码矩阵计算
MMSE_matrix=zeros(M,N,user_num,L);
for k=1:user_num
    for s=1:L
        H_UL_ks=Channel_mul((1+(k-1)*N):k*N,:,s,s)';
        R_F=zeros(M,M);
        for i=1:user_num
        for l=1:L
            H_UL_il=Channel_mul((1+(i-1)*N):i*N,:,s,l)';
            R_F=R_F+power*H_UL_il*H_UL_il';       
        end    
        end
        R_F=R_F+var_noise*eye(M);
        MMSE_matrix(:,:,k,s)=pinv(R_F)*sqrt(power)*H_UL_ks;
        
    end
    
end

%%计算SE
MMSE_rates=zeros(1,user_num);
for k=1:user_num
     SE_k=zeros(1,L);
    for s=1:L
        w_k_s=MMSE_matrix(:,:,k,s)/norm(MMSE_matrix(:,:,k,s),'fro');
        H_ss=Channel_mul((1+(k-1)*N):k*N,:,s,s);
        B_kk=sqrt(power)*H_ss*w_k_s;
        
        inter=zeros(N,N);
        
        for i=1:user_num
            if i~=k
                w_i_s=MMSE_matrix(:,:,i,s)/norm(MMSE_matrix(:,:,i,s),'fro');
                B_ki_s=sqrt(power)*H_ss*w_i_s;
                inter=inter+B_ki_s*B_ki_s';
            end
        end
        
        for i=1:user_num
              w_i_s=MMSE_matrix(:,:,i,s)/norm(MMSE_matrix(:,:,i,s),'fro');
            for l=1:L
                if l~=s
                   H_sl=Channel_mul((1+(k-1)*N):k*N,:,s,l);
                   B_ki_l=sqrt(power)*H_sl*w_i_s;
                    inter=inter+B_ki_l*B_ki_l';
                end      
            end    
        end
        
        inter=inter+var_noise*eye(N);
        SE_k(s)=log2(det(eye(N)+B_kk'*inv(inter)*B_kk));      
    end
    MMSE_rates(k)=sum(SE_k);
end

MMSE=sum( MMSE_rates)/user_num;
end