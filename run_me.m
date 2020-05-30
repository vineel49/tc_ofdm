% Turbo coded (2x2) MIMO-OFDM

clear all
close all
clc
num_bit = 1024; % number of data bits, overall rate 1/2
FFT_len = 1024; % length of the fft at each transmit antenna
SNR_dB = 0; % SNR per bit in dB (in logarithmic scale)
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
ofdm_intr_map = randperm(num_bit);
ofdm_deintr_map = deintrlv((1:num_bit),ofdm_intr_map);

%-------------Defining SNR per bit - overall rate is 1/2 ------------------
SNR = 10^(0.1*SNR_dB); % SNR is average SNR per bit (linear scale)
NOISE_VAR_1D = 2*fade_var_1D*chan_len*2*2*2/(2*FFT_len*SNR*1*2); % 1D noise variance

% QPSK Symbols
QPSK_SYM = zeros(4,num_bit);
QPSK_SYM(1,:) = (1+1i)*ones(1,num_bit);
QPSK_SYM(2,:) = (1-1i)*ones(1,num_bit);
QPSK_SYM(3,:) = (-1+1i)*ones(1,num_bit);
QPSK_SYM(4,:) = (-1-1i)*ones(1,num_bit);

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

% interleaving to decorrelate the channel
mod_sig1 = mod_sig1(ofdm_intr_map);
mod_sig2 = mod_sig2(ofdm_intr_map);

% IFFT operation
F_trans_sig_no_CP1 = ifft(mod_sig1);
F_trans_sig_no_CP2 = ifft(mod_sig2);

% inserting cyclic prefix
F_trans_sig1 = [F_trans_sig_no_CP1(end-cp_len+1:end) F_trans_sig_no_CP1];
F_trans_sig2 = [F_trans_sig_no_CP2(end-cp_len+1:end) F_trans_sig_no_CP2];
%----------------  CHANNEL -----------------------------------------
% Rayleigh channel
% channel for 1st tx, 1st rx 
fade_chan11 = sqrt(fade_var_1D)*randn(1,chan_len) + 1i*sqrt(fade_var_1D)*randn(1,chan_len);     
F_fade_chan11 = fft(fade_chan11,FFT_len);
deintr_F_fade_chan11 = F_fade_chan11(ofdm_deintr_map);


% channel for 2nd tx, 1st rx 
fade_chan21 = sqrt(fade_var_1D)*randn(1,chan_len) + 1i*sqrt(fade_var_1D)*randn(1,chan_len);     
F_fade_chan21 = fft(fade_chan21,FFT_len);
deintr_F_fade_chan21 = F_fade_chan21(ofdm_deintr_map);

% channel for 1st tx, 2nd rx 
fade_chan12 = sqrt(fade_var_1D)*randn(1,chan_len) + 1i*sqrt(fade_var_1D)*randn(1,chan_len);     
F_fade_chan12 = fft(fade_chan12,FFT_len);
deintr_F_fade_chan12 = F_fade_chan12(ofdm_deintr_map);

% channel for 2nd tx, 2nd rx 
fade_chan22 = sqrt(fade_var_1D)*randn(1,chan_len) + 1i*sqrt(fade_var_1D)*randn(1,chan_len);     
F_fade_chan22 = fft(fade_chan22,FFT_len);
deintr_F_fade_chan22 = F_fade_chan22(ofdm_deintr_map);

% noise at rx 1
noise1 = sqrt(NOISE_VAR_1D)*randn(1,FFT_len+cp_len + chan_len-1) + 1i*sqrt(NOISE_VAR_1D)*randn(1,FFT_len+cp_len + chan_len-1);    

% noise at rx 2
noise2 = sqrt(NOISE_VAR_1D)*randn(1,FFT_len+cp_len + chan_len-1) + 1i*sqrt(NOISE_VAR_1D)*randn(1,FFT_len+cp_len + chan_len-1);    

% channel output
% rx 1
Chan_Op1 = conv(F_trans_sig1,fade_chan11)+conv(F_trans_sig2,fade_chan21) + noise1;  
% rx 2
Chan_Op2 = conv(F_trans_sig1,fade_chan12)+conv(F_trans_sig2,fade_chan22) + noise2; 
%-----------------   RECEIVER  --------------------------------------------
% rx 1
Chan_Op1(1:cp_len) = [];
Rec_Sig1 = Chan_Op1(1:FFT_len);

% at rx 2
Chan_Op2(1:cp_len) = [];
Rec_Sig2 = Chan_Op2(1:FFT_len);

% fft operation at rx 1
F_rec_sig1 = fft(Rec_Sig1);
% deinterleaving
deintr_F_rec_sig1 = F_rec_sig1(ofdm_deintr_map);

% fft operation at rx 2
F_rec_sig2 = fft(Rec_Sig2);
deintr_F_rec_sig2 = F_rec_sig2(ofdm_deintr_map);

Dist = zeros(4,4,FFT_len);
% Branch metrices for the SOVA
  for tx1_sym=1:4
   for tx2_sym=1:4
    Dist(tx1_sym,tx2_sym,:)=abs(deintr_F_rec_sig1...
                           -QPSK_SYM(tx1_sym,:).*deintr_F_fade_chan11...
                           -QPSK_SYM(tx2_sym,:).*deintr_F_fade_chan21).^2 ...
                           +abs(deintr_F_rec_sig2...
                           -QPSK_SYM(tx1_sym,:).*deintr_F_fade_chan12...
                           -QPSK_SYM(tx2_sym,:).*deintr_F_fade_chan22).^2;
   end
  end

 [M1,~]=min(Dist,[],2);
 [M2,~]=min(Dist,[],1); 
  
 branch_metric1 = reshape(M1,4,FFT_len);% branch metrices for decoder 1
 branch_metric2 = reshape(M2,4,FFT_len);% branch metrices for decoder 2


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
soft_output = SOVA_END(soft_output(intr_map),num_bit,branch_metric2); %8

% hard decision is taken on the a posteriori probabilities
soft_output  = soft_output(deintr_map);
dec_a = soft_output<0;
C_Ber = C_Ber + nnz(a-dec_a);
end
% Bit error rate
BER =  C_Ber/(num_bit*num_frames)

toc()


