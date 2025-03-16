function[RecLength,ChannelSigma,RecResVector,TraResVector,nr_act,ns_act]=channel_generate(laisi_k,UserSpacing,RecSpacing,TraSpacing,lambda,RecLength,TraLength,UsersNum,double)

%% Generate the channel matrix
% Obtain the sampled eclipse and their coordinates given a surface
[RecSamplePoint,WD_RecX_vec,WD_RecY_vec]=WD_eclipse_sample(lambda,RecLength);
[TraSamplePoint,WD_TraX_vec,WD_TraY_vec]=WD_eclipse_sample(lambda,TraLength);

[RecVarianceVec,WD_RecX_vec,WD_RecY_vec,nr_act]=WaveDomainChannel(RecSamplePoint,WD_RecX_vec,WD_RecY_vec,lambda,RecLength);
[TraVarianceVec,WD_TraX_vec,WD_TraY_vec,ns_act]=WaveDomainChannel(TraSamplePoint,WD_TraX_vec,WD_TraY_vec,lambda,TraLength);

% average
if double~=1
    tt1=mean(TraVarianceVec);
    TraVarianceVec=tt1*ones(ns_act,1);
end
ChannelSigma_single=RecVarianceVec*TraVarianceVec';
if UsersNum>0
    ChannelSigma=repmat(ChannelSigma_single,UsersNum,1);
    RecVarianceVec=repmat(RecVarianceVec,UsersNum,1);
else
    ChannelSigma=ChannelSigma_single;
end

[TraResVector,TraNumNs]=ResponseGenerate(TraSpacing,TraSpacing,lambda,TraLength,WD_TraX_vec,WD_TraY_vec,'Transmit'); 

% Multi-user Case
RecResVector=[];
for user=1:UsersNum 
    RecLength(user)=RecLength(1);
    RecLength(user).distance(1) =(user-1)*UserSpacing+RecLength(1).distance(1);
    tempUserInfo=RecLength(user);
    distances=norm(RecLength(user).distance-TraLength.distance,2);
         [channelGaindB,~] = functionChannelgain(distances);
        channelGain=sqrt(1/(1+laisi_k))*sqrt(db2pow(channelGaindB)); 
    [RecResVectorPart,RecNumNr]=ResponseGenerate(RecSpacing,RecSpacing,lambda,tempUserInfo,WD_RecX_vec,WD_RecY_vec,'Receive');
    RecResVector=blkdiag(RecResVector,RecResVectorPart);
end
 RecResVector=channelGain*RecResVector;
%% Pass Channel 
% [Phi]=PhaseGenerate(TraNumNs);  % Phase matrix 
Phi=eye(TraNumNs); % Phi is an identity matrix
TransSig=randn(UsersNum*RecNumNr,1)+1i*randn(UsersNum*RecNumNr,1); % transmitted signal
%%%生成NLOS信道的采样点

end