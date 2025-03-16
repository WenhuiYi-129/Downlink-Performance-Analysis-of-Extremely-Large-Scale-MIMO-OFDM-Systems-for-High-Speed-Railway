function [MRT_rates,MRT]=BFSchemes_MRT(Channel_mul,var_noise,power,user_num,L,N,M)
 
%% %MRT for Perfect channels
% MRT_rates=zeros(user_num,L);
%%预编码矩阵计算
MRT_matrix=zeros(M,N,user_num,L);
for k=1:user_num
    for s=1:L
        H_UL_ks=Channel_mul((1+(k-1)*N):k*N,:,s,s)';
        MRT_matrix(:,:,k,s)=H_UL_ks;   
    end   
end
%%计算SE
MRT_rates=zeros(1,user_num);
for k=1:user_num
     SE_k=zeros(1,L);
    for s=1:L
        w_k_s=MRT_matrix(:,:,k,s)/norm(MRT_matrix(:,:,k,s),'fro');
        H_ss=Channel_mul((1+(k-1)*N):k*N,:,s,s);
        B_kk=sqrt(power)*H_ss*w_k_s;
        
        inter=zeros(N,N);
        
        for i=1:user_num
            if i~=k
                w_i_s=MRT_matrix(:,:,i,s)/norm(MRT_matrix(:,:,i,s),'fro');
                B_ki_s=sqrt(power)*H_ss*w_i_s;
                inter=inter+B_ki_s*B_ki_s';
            end
        end
        
        for i=1:user_num
              w_i_s=MRT_matrix(:,:,i,s)/norm(MRT_matrix(:,:,i,s),'fro');
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
    MRT_rates(k)=sum(SE_k);
end

MRT=sum( MRT_rates)/user_num;
end