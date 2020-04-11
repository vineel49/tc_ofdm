% ideal coherent detection - sova 
clear all
close all
clc
num_bit = 512; % number of data bits, overall rate 1/2
FFT_len = 1024; % length of the fft
SNR_dB = 9; % SNR per bit in dB (in logarithmic scale)
fade_var_1D= 0.5; % 1D fade variance
num_frames = 10^1; % number of frames
chan_len = 10; % number of channel taps
cp_len = chan_len-1; % length of the cyclic prefix

%------------- GENERATOR polynomial of the encoder-------------------------
gen_poly = ldiv2([1 0 1],[1 1 1],num_bit); % using long division method

%------------- Interleaver mapping of the turbo code ----------------------
intr_map = randperm(num_bit);
deintr_map = deintrlv((1:num_bit),intr_map);

%------------- ofdm interleaver--------------------------------------------
ofdm_intr_map = randperm(2*num_bit);

%-------------Defining SNR per bit - overall rate is 1/2 ------------------
SNR = 10^(0.1*SNR_dB); % SNR is average SNR per bit (linear scale)
NOISE_VAR_1D = 2*2*2*fade_var_1D*chan_len/(2*FFT_len*SNR); % 1D noise variance

C_Ber = 0;
tic()
for frame_cnt = 1:num_frames
%-----------------     TRANSMITTER    -------------------------------------
% source bit generation
a= randi([0 1],1,num_bit);

% Turbo encoder
% encoder 1
b1 = zeros(1,2*num_bit); % encoder 1 output initialization
b1(1:2:end) = a; % systematic bit
temp1 = mod(conv(gen_poly,a),2); 
b1(2:2:end) = temp1(1:num_bit); % parity bit
% encoder 2
b2 = zeros(1,2*num_bit); % encoder 2 output initialization
b2(1:2:end) = a(intr_map); % systematic bit
temp2 = mod(conv(gen_poly,b2(1:2:end)),2); 
b2(2:2:end) = temp2(1:num_bit); % parity bit

% QPSK modulation (according to set partitioning principles)
mod_sig1 = 1-2*b1(1:2:end) + 1i*(1-2*b1(2:2:end));
mod_sig2 = 1-2*b2(1:2:end) + 1i*(1-2*b2(2:2:end));
mod_sig = [mod_sig1 mod_sig2];

% interleaving to decorrelate the channel
mod_sig = mod_sig(ofdm_intr_map);

% IFFT operation
F_trans_sig_no_CP = ifft(mod_sig);

% inserting cyclic prefix
F_trans_sig = [F_trans_sig_no_CP(end-cp_len+1:end) F_trans_sig_no_CP];
%----------------  CHANNEL -----------------------------------------
% Rayleigh channel
fade_chan = sqrt(fade_var_1D)*randn(1,chan_len) + 1i*sqrt(fade_var_1D)*randn(1,chan_len);     
F_fade_chan = fft(fade_chan,FFT_len);
noise = sqrt(NOISE_VAR_1D)*randn(1,FFT_len+cp_len + chan_len-1) + 1i*sqrt(NOISE_VAR_1D)*randn(1,FFT_len+cp_len + chan_len-1);    
Chan_Op = conv(F_trans_sig,fade_chan) + noise;    
%-----------------   RECEIVER  --------------------------------------------
Chan_Op(1:cp_len) = [];
Rec_Sig = Chan_Op(1:FFT_len);

% fft operation
F_rec_sig = fft(Rec_Sig);


% Branch metrices for the SOVA
QPSK_SYM = zeros(4,2*num_bit);
QPSK_SYM(1,:) = (1+1i)*ones(1,2*num_bit);
QPSK_SYM(2,:) = (1-1i)*ones(1,2*num_bit);
QPSK_SYM(3,:) = (-1+1i)*ones(1,2*num_bit);
QPSK_SYM(4,:) = (-1-1i)*ones(1,2*num_bit);

branch_metric = zeros(4,2*num_bit);
 branch_metric(1,ofdm_intr_map(:))=abs(F_rec_sig-QPSK_SYM(1,:).*F_fade_chan).^2;
 branch_metric(2,ofdm_intr_map(:))=abs(F_rec_sig-QPSK_SYM(2,:).*F_fade_chan).^2;
 branch_metric(3,ofdm_intr_map(:))=abs(F_rec_sig-QPSK_SYM(3,:).*F_fade_chan).^2;
 branch_metric(4,ofdm_intr_map(:))=abs(F_rec_sig-QPSK_SYM(4,:).*F_fade_chan).^2;
 
branch_metric1 = branch_metric(:,1:num_bit); % branch metrices for decoder 1
branch_metric2 = branch_metric(:,num_bit+1:end); % branch metrices for decoder 2
 
% a priori probabilities (LLR) - initialization
apr_LLR = zeros(1,num_bit); % for first iteration

% iterative decoding
soft_output = SOVA(apr_LLR,num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %1

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %2

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %3

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %4

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %5

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %6

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %7

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %8

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA(soft_output(intr_map),num_bit,branch_metric2); %9

soft_output = SOVA(soft_output(deintr_map),num_bit,branch_metric1);
soft_output = SOVA_END(soft_output(intr_map),num_bit,branch_metric2); %10

% hard decision is taken on the a posteriori probabilities
soft_output  = soft_output(deintr_map);
dec_a = soft_output<0;
C_Ber = C_Ber + nnz(a-dec_a);
end
% Bit error rate
BER =  C_Ber/(num_bit*num_frames)

toc()


