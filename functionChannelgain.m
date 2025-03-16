function [channelGaindB,ricianFactor] = functionChannelgain(distances)

%Model parameters

%Standard deviation of shadow fading in dB

sigma_sf_LOS=4;   %for LOS

%Prepare to store channel gain (in dB), Rician Factor \kappa and
%probabaility of LoS for each UE
% channelGaindB=zeros(K,L,L);
% probLOS=zeros(K,L,L);

 ricianFactor=db2pow(13-0.03*distances);
        %Compute average channel gain using the large-scale fading
        %model based on 3GPP model while neglecting the shadow fading for each UE

   channelGaindB= -30.18-26*log10(distances);
    %Go through all UEs in cell l and generate shadow fading realizations        
            shadowing = sigma_sf_LOS*randn(1);
            channelGainShadowing = channelGaindB+ shadowing;
        %Store average channel gains with shadowing fading
        channelGaindB= channelGainShadowing;

    
    

