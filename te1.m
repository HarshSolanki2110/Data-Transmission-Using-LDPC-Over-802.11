clc
clear
close all
maxnumiter = 10;
% snr = linspace(1,7,10);
snr = 1:10;
numframes = 10000;
P = Proto_mat_Gen();
blockSize =96;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);
figure(WindowState="maximized")
spy(pcmatrix)
grid("on")
grid("minor")
cfgLDPCEnc = ldpcEncoderConfig(pcmatrix);
cfg_bp_LDPCDec = ldpcDecoderConfig(pcmatrix,'bp');
cfg_lbp_LDPCDec = ldpcDecoderConfig(pcmatrix,"layered-bp");
cfg_nmbp_LDPCDec = ldpcDecoderConfig(pcmatrix,"norm-min-sum");
cfg_onmbp_LDPCDec = ldpcDecoderConfig(pcmatrix,"offset-min-sum");
bpskmod = comm.BPSKModulator;
% bpskmod.PhaseOffset = pi/16;
bpskdemod = comm.BPSKDemodulator;
bpskdemod.DecisionMethod = 'Log-likelihood ratio';
EbN0 = 10.^(snr/10);
BER_NoCoding = 1/2.*erfc(sqrt(EbN0));
BER_bp = zeros(1,size(snr,2));
BER_lbp = zeros(1,size(snr,2));
BER_nmbp = zeros(1,size(snr,2));
BER_onmbp = zeros(1,size(snr,2));
FER_bp = 0;
FER_lbp = 0;
FER_nmbp = 0;
FER_onmbp = 0;
for ii = 1:length(snr)
    errNoCoding = 0;
    err_bp = 0;
    err_lbp = 0;
    err_nmbp = 0;
    err_onmbp = 0;
    F_err_bp = 0;
    F_err_lbp = 0;
    F_err_nmbp = 0;
    F_err_onmbp = 0;
    for counter = 1:numframes
        data = randi([0 1],cfgLDPCEnc.NumInformationBits,1);
        % Transmit and receive with LDPC coding
        encodedData = ldpcEncode(data,cfgLDPCEnc);
%         modSignal = bpskmod(encodedData);
        modSignal = BPSK_Calcu(EbN0(ii),encodedData);
        receivedSignal = awgn(modSignal,EbN0(ii));
        demodSignal = bpskdemod(receivedSignal);
        received_bp_Bits = ldpcDecode(demodSignal,cfg_bp_LDPCDec,maxnumiter,"DecisionType","hard");
        received_lbp_Bits = ldpcDecode(demodSignal,cfg_lbp_LDPCDec,maxnumiter,"DecisionType","hard");
        received_nmbp_Bits = ldpcDecode(demodSignal,cfg_nmbp_LDPCDec,maxnumiter,"DecisionType","hard");
        received_onmbp_Bits = ldpcDecode(demodSignal,cfg_onmbp_LDPCDec,maxnumiter,"DecisionType","hard");
        err_bp_B = sum(data' ~= received_bp_Bits',2);
        err_lbp_B = sum(data' ~= received_lbp_Bits',2);
        err_nmbp_B = sum(data' ~= received_nmbp_Bits',2);
        err_onmbp_B = sum(data' ~= received_onmbp_Bits',2);
        if err_bp_B ==0
            F_err_bp = F_err_bp + 1;
        end
        if err_lbp_B ==0
            F_err_lbp = F_err_lbp + 1;
        end
        if err_nmbp_B ==0
            F_err_nmbp = F_err_nmbp + 1;
        end
        if err_onmbp_B ==0
            F_err_onmbp = F_err_onmbp + 1;
        end
        err_bp = err_bp + err_bp_B;
        err_lbp = err_lbp + err_lbp_B;
        err_nmbp = err_nmbp + err_nmbp_B;
        err_onmbp = err_onmbp + err_onmbp_B;
    end
%     BER_NoCoding(:,ii) = errNoCoding/(cfg_bp_LDPCDec.NumInformationBits*numframes*size(snr,2));
    BER_bp(:,ii) = err_bp/(cfg_bp_LDPCDec.NumInformationBits*numframes*size(snr,2));
    BER_lbp(:,ii) = err_lbp/(cfg_bp_LDPCDec.NumInformationBits*numframes*size(snr,2));
    BER_nmbp(:,ii) = err_nmbp/(cfg_bp_LDPCDec.NumInformationBits*numframes*size(snr,2));
    BER_onmbp(:,ii) = err_onmbp/(cfg_bp_LDPCDec.NumInformationBits*numframes*size(snr,2));
    FER_bp(:,ii) = F_err_bp/numframes;
    FER_lbp(:,ii) = F_err_lbp/numframes;
    FER_nmbp(:,ii) = F_err_nmbp/numframes;
    FER_onmbp(:,ii) = F_err_onmbp/numframes;
end
figure(WindowState="maximized")
semilogy(snr,BER_NoCoding,"Marker","*",LineWidth=1);
hold on
semilogy(snr,BER_bp,"Marker","^",LineWidth=1);
semilogy(snr,BER_lbp,"Marker","diamond",LineWidth=1);
semilogy(snr,BER_nmbp,"Marker","o",LineWidth=1);
semilogy(snr,BER_onmbp,"Marker","hexagram",LineWidth=1);
grid("on")
grid("minor")
% xlim([-1 10]);
% ylim([10e-6 10e-4]);
xlabel ('Signal to Noise Ratio in dB',FontSize=12,Color=[1,0,1])
ylabel('Bit Error Rate',FontSize=12,Color=[1,0,1])
legend("No Codeing","BP Decoder","Linear BP","Normalized Min-Sum BP","Offset Min-Sum BP","Location","best")
savefig('BER_All_1_2_1152.fig')
saveas(gcf,'BER_All_1_2_1152.png','png')
saveas(gcf,'BER_All_1_2_1152.jpeg','jpeg')


figure(WindowState="maximized")
semilogy(snr,FER_bp,"Marker","^",LineWidth=1);
hold on
semilogy(snr,FER_lbp,"Marker","diamond",LineWidth=1);
semilogy(snr,FER_nmbp,"Marker","o",LineWidth=1);
semilogy(snr,FER_onmbp,"Marker","hexagram",LineWidth=1);
grid("on")
grid("minor")
% xlim([-1 10]);
% ylim([10e-6 10e-4]);
xlabel ('Signal to Noise Ratio in dB',FontSize=12,Color=[1,0,1])
ylabel('Frame Error Rate',FontSize=12,Color=[1,0,1])
legend("FER BP Decoder","FER Linear BP","FER Normalized Min-Sum BP","FER Offset Min-Sum BP","Location","best")
BER_All = [BER_NoCoding;BER_bp;BER_lbp;BER_nmbp;BER_onmbp];
FER_All = [FER_bp;FER_lbp;FER_nmbp;FER_onmbp];
savefig('FER_All_1_2_1152.fig')
saveas(gcf,'FER_All_1_2_1152.png','png')
saveas(gcf,'FER_All_1_2_1152.jpeg','jpeg')