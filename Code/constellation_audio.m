clc, clearvars, close all;
%A/D converter
[audio, fs] = audioread('project.wav'); % Load your audio file here
audio_normalized = int16(audio * 32767); % Normalize
audio_binary = dec2bin(typecast(audio_normalized(:), 'uint16'), 16); % Convert to binary
binary_vector = audio_binary(:)';
numeric_vector = double(binary_vector) - 48; % Reshape to vector for transmission
%numeric_vector=numeric_vector(1:1300000);
% numeric_vector = [0,1,1,0,0,0,1,1];
% Encoder
j = 1;
mapped_symbols = zeros(1, (length(numeric_vector)/2));
for i = 1:2:length(numeric_vector) - 1
    bit1 = numeric_vector(i);
    bit2 = numeric_vector(i+1);
    mapped_symbols(j) = qpskmap(bit1, bit2);
    j = j+1;
end
mapped_symbols_real=real(mapped_symbols);
mapped_symbols_imag=imag(mapped_symbols);
% figure;
% plot(real(mapped_symbols), imag(mapped_symbols), 'bo');
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('QPSK Constellation');
% ylim([-3 3]);
% xlim([-3 3]);
% figure;
% stem(mapped_symbols_real(1:200), 'bo');
% title('mapped sym real');
% figure;
% stem(mapped_symbols_imag(1:200), 'bo');
% title('mapped sym imag');
%Line coding
oversampling_factor = 9;
m = oversampling_factor;
m2=10;
a = 0.5;
length_ = 20;
[transmit_filter_rc, ~] = raised_cosine(a, m, 9);
nsymbols_upsampled_r = 1 + (length(mapped_symbols_real) - 1) * m;
symbols_upsampled_r = zeros(nsymbols_upsampled_r, 1);
symbols_upsampled_r(1:m:nsymbols_upsampled_r) = mapped_symbols_real;
nsymbols_upsampled_imag = 1 + (length(mapped_symbols_imag) - 1) * m;
symbols_upsampled_imag = zeros(nsymbols_upsampled_imag, 1);
symbols_upsampled_imag(1:m:nsymbols_upsampled_imag) = mapped_symbols_imag;
nsymbols_upsampled_r_rect = 1 + (length(mapped_symbols_real) - 1) * m2;
symbols_upsampled_r_rect = zeros(nsymbols_upsampled_r_rect, 1);
symbols_upsampled_r_rect(1:m2:nsymbols_upsampled_r_rect) = mapped_symbols_real;
nsymbols_upsampled_imag_rect = 1 + (length(mapped_symbols_imag) - 1) * m2;
symbols_upsampled_imag_rect = zeros(nsymbols_upsampled_imag_rect, 1);
symbols_upsampled_imag_rect(1:m2:nsymbols_upsampled_imag_rect) = mapped_symbols_imag;
% figure;
% subplot(2,1,1);
% plot(symbols_upsampled_imag(1:1000), 'bo');
% subplot(2,1,2);
% plot(symbols_upsampled_r(1:1000));
% title("upsampled data");
rect_pulse = [ones(1, floor(length_/2)), zeros(1, floor(length_/2))];
transmit_filter_rect = rect_pulse;
% figure;
% stem(rect_pulse);
tx_output_rc_r = conv(symbols_upsampled_r, transmit_filter_rc, 'same');
tx_output_rc_imag = conv(symbols_upsampled_imag, transmit_filter_rc, 'same');
tx_output_rect_r = conv(symbols_upsampled_r_rect, rect_pulse, 'same');
tx_output_rect_imag = conv(symbols_upsampled_imag_rect, rect_pulse, 'same');
% tx_output_rect = symbols_upsampled;
% figure;
% plot(real(tx_output_rc), imag(tx_output_rc), 'bo');
% title("tx output rc");
% 
% figure;
% stem(tx_output_rect_r(1:1000), 'ro');
% title("upsampled rect and conv for rect real");
% figure;
% stem(tx_output_rect_imag(1:1000), 'ro');
% title("upsampled rect and conv for rect imag");
% figure;
% stem(tx_output_rc_r(1:1000), 'ro');
% title("upsampled rect and conv for rc real");
% figure;
% stem(tx_output_rc_imag(1:1000), 'ro');
% title("upsampled rect and conv for rc imag");
% Modulation_
fc = 1e6;
%t = (0:length(tx_output_rc_r)-1) / fs;
t=0:1/(6*fc):(length(tx_output_rc_r)-1)/(6*fc);
% t_rect = (0:length(tx_output_rect_imag)-1) / fs;
t_rect = 0:1/(6*fc):(length(tx_output_rect_imag)-1)/(6*fc);
% % Limit the number of values for modulation
% limit = 1000;
% if length(tx_output_rc) > limit
%     tx_output_rc = tx_output_rc(100002:100100);
%     t = t(1:99);
% end
% 
% if length(tx_output_rect) > limit
%     tx_output_rect = tx_output_rect(100002:100100);
%     t_rect = t_rect(1:99);
% end
% Modulate tx_output_rc
tx_output_rc_half1 = tx_output_rc_r;
tx_output_rc_half2 = tx_output_rc_imag;
modulated_half1_rc = tx_output_rc_half1 .* cos(2*pi*fc*t)';
modulated_half2_rc = tx_output_rc_half2 .* sin(2*pi*fc*t)';
modulated_signal_rc = modulated_half1_rc + modulated_half2_rc;
% Modulate tx_output_rect
tx_output_rect_half1 = tx_output_rect_r;
tx_output_rect_half2 = tx_output_rect_imag;
modulated_half1_rect = (tx_output_rect_half1) .* cos(2*pi*fc*t_rect)';
modulated_half2_rect = tx_output_rect_half2 .* sin(2*pi*fc*t_rect)';
modulated_signal_rect = modulated_half1_rect + modulated_half2_rect;
% modulated_signal_rect = awgn(modulated_signal_rect,8,"measured");
% figure;
% subplot(2,1,1);
% plot(modulated_half1_rect(630:1000));
% title("tx output rect real")
% subplot(2,1,2);
% plot(modulated_half2_rect(630:1000));
% title("tx output rect imag")
% figure;
% plot(modulated_signal_rect(630:1000));
% title("modulated signal");
% figure;
% subplot(2,1,1);
% plot(modulated_half1_rc(630:1000));
% title("tx output rc real")
% subplot(2,1,2);
% plot(modulated_half2_rc(630:1000));
% title("tx output rc imag")
% figure;
% plot(modulated_signal_rc(630:1000));
% title("modulated signal for rc");
% % with
% noise..................................................................................
a=0.1;
b=1;
tb = log2(4)*1/fs;  
%t_temp = 1:1:(length(modulated_signal_rect));
%t_temp_rc = 1:1:(length(modulated_signal_rc));
h_rect=zeros(length(modulated_signal_rect),1);
h_rc=zeros(length(modulated_signal_rc),1);
for t_temp=1:1:(length(modulated_signal_rect))
    if(t_temp-b*tb>0)
    h_rect = a*modulated_signal_rect(t_temp) + (1 - a)*modulated_signal_rect(round(t_temp - b*tb));
    end
end
for t_temp_rc=1:1:(length(modulated_signal_rc))
    if(t_temp_rc-b*tb>0)
    h_rc = a*modulated_signal_rc(t_temp_rc) + (1 - a)*modulated_signal_rc(round(t_temp_rc - b*tb));
    end
end
n_awgn_rect = awgn(modulated_signal_rect,8,"measured");
n_awgn_rc = awgn(modulated_signal_rc,8,"measured");
% rx_memoryless_rc=modulated_signal_rc+n_awgn_rc;
% rx_memoryless_retc=modulated_signal_rect+n_awgn_rect;
% 
 rx_memory_rect = h_rect+ n_awgn_rect;
 rx_memory_rc = h_rc+ n_awgn_rc;
% plot(n_awgn_rc);
% title("modulated rc with noise");
% plot(n_awgn_rect);
% title("modulated rect with noise");
%demodulator
%.........................................................................................
demodu=zeros(length(numeric_vector),1);
demodu_rc=zeros(length(numeric_vector),1);
Rx_in=0;
Rx_qd=0;
Rx_in_rc=0;
Rx_qd_rc=0;
% for i=1:1:ceil(length(modulated_signal_rect)/2)
    % start_index = (i - 1)*2  + 1;
    % end_index = min(i*2 , length(modulated_signal_rect));
    % t_rect_demod = start_index:2/fs:end_index;
    x_in=n_awgn_rect.*cos(2*pi*fc*t_rect)';
    x_in_rc=n_awgn_rc.*cos(2*pi*fc*t)';
    % %in_intg=(trapz(t,x_in))*(2/T);
    % in_intg=(trapz(1,x_in));
    x_qd=n_awgn_rect.*sin(2*pi*fc*t_rect)';
    x_qd_rc=n_awgn_rc.*sin(2*pi*fc*t)';
    %qd_intg=(trapz(t,x_qd))*(2/T);    
    % for k = 1:length(x_in)
    %   qd_intg=(trapz(2/fs,x_qd(i)));
    %   in_intg=(trapz(2/fs,x_in(i)));
    % end
    % qd_intg=(trapz(1,x_qd));
    line_decoder_rect_qd = conv(x_qd, transmit_filter_rect, 'same');
     downsample_rect_qd=downsample(line_decoder_rect_qd,10);
%
    line_decoder_rect_in = conv(x_in, transmit_filter_rect, 'same');
     downsample_rect_in=downsample(line_decoder_rect_in,10);
     line_decoder_rc_qd = conv(x_qd_rc, transmit_filter_rc, 'same');
     downsample_rc_qd=downsample(line_decoder_rc_qd,9);
%
    line_decoder_rc_in = conv(x_in_rc, transmit_filter_rc, 'same');
     downsample_rc_in=downsample(line_decoder_rc_in,9);
%
j = 1;
for i = 1:length(numeric_vector)/2
    if(downsample_rect_in(i)<0 && downsample_rect_qd(i)<0)
        Rx_in = 1;
        Rx_qd = 1;
    elseif(downsample_rect_in(i)<0 && downsample_rect_qd(i)>0)
        Rx_in = 1;
        Rx_qd = 0;
    elseif(downsample_rect_in(i)>0 && downsample_rect_qd(i)>0)
        Rx_in = 0;
        Rx_qd = 0;
    elseif(downsample_rect_in(i)>0 && downsample_rect_qd(i)<0)
        Rx_in = 0;
        Rx_qd = 1;
    end
    demodu(j)=Rx_in;
    demodu(j+1) = Rx_qd;
       j=j+2;
end
k = 1;
for z = 1:length(numeric_vector)/2
    if(downsample_rc_in(z)<0 && downsample_rc_qd(z)<0)
        Rx_in_rc = 1;
        Rx_qd_rc = 1;
    elseif(downsample_rc_in(z)<0 && downsample_rc_qd(z)>0)
        Rx_in_rc = 1;
        Rx_qd_rc = 0;
    elseif(downsample_rc_in(z)>0 && downsample_rc_qd(z)>0)
        Rx_in_rc = 0;
        Rx_qd_rc = 0;
    elseif(downsample_rc_in(z)>0 && downsample_rc_qd(z)<0)
        Rx_in_rc = 0;
        Rx_qd_rc = 1;
    end
    demodu_rc(k)=Rx_in_rc;
    demodu_rc(k+1) = Rx_qd_rc;
       k=k+2;
end
% figure();
% subplot(2,1,1);
% plot(numeric_vector(1:9000));
% title("orig for rect");
% subplot(2,1,2);
% plot(demodu(1:9000));
% title("demodulated for rect");
% figure();
% subplot(2,1,1);
% plot(numeric_vector(1:9000));
% title("orig for rc");
% subplot(2,1,2);
% stem(demodu_rc(1:9000));
% title("demodulated for rc");
demodu_m=zeros(length(numeric_vector),1);
demodu_m_rc=zeros(length(numeric_vector),1);
Rx_in_m=0;
Rx_qd_m=0;
Rx_in_m_rc=0;
Rx_qd_m_rc=0;
% for i=1:1:ceil(length(modulated_signal_rect)/2)
    % start_index = (i - 1)*2  + 1;
    % end_index = min(i*2 , length(modulated_signal_rect));
    % t_rect_demod = start_index:2/fs:end_index;
    x_in_m=rx_memory_rect.*cos(2*pi*fc*t_rect)';
    x_in_rc_m=rx_memory_rc.*cos(2*pi*fc*t)';
    % %in_intg=(trapz(t,x_in_m))*(2/T);
    % in_intg=(trapz(1,x_in_m));
    x_qd_m=rx_memory_rect.*sin(2*pi*fc*t_rect)';
    x_qd_rc_m=rx_memory_rc.*sin(2*pi*fc*t)';
    %qd_intg=(trapz(t,x_qd_m))*(2/T);    
    % for k = 1:length(x_in_m)
    %   qd_intg=(trapz(2/fs,x_qd_m(i)));
    %   in_intg=(trapz(2/fs,x_in_m(i)));
    % end
    % qd_intg=(trapz(1,x_qd_m));
    line_decoder_rect_qd_m = conv(x_qd_m, transmit_filter_rect, 'same');
     downsample_rect_qd_m=downsample(line_decoder_rect_qd_m,10);
%
    line_decoder_rect_in_m = conv(x_in_m, transmit_filter_rect, 'same');
     downsample_rect_in_m=downsample(line_decoder_rect_in_m,10);
     line_decoder_rc_qd_m = conv(x_qd_rc_m, transmit_filter_rc, 'same');
     downsample_rc_qd_m=downsample(line_decoder_rc_qd_m,9);
%
    line_decoder_rc_in_m = conv(x_in_rc_m, transmit_filter_rc, 'same');
     downsample_rc_in_m=downsample(line_decoder_rc_in_m,9);
%
v = 1;
for i = 1:length(numeric_vector)/2
    if(downsample_rect_in_m(i)<0 && downsample_rect_qd_m(i)<0)
        Rx_in_m = 1;
        Rx_qd_m = 1;
    elseif(downsample_rect_in_m(i)<0 && downsample_rect_qd_m(i)>0)
        Rx_in_m = 1;
        Rx_qd_m = 0;
    elseif(downsample_rect_in_m(i)>0 && downsample_rect_qd_m(i)>0)
        Rx_in_m = 0;
        Rx_qd_m = 0;
    elseif(downsample_rect_in_m(i)>0 && downsample_rect_qd_m(i)<0)
        Rx_in_m = 0;
        Rx_qd_m = 1;
    end
    demodu_m(v)=Rx_in_m;
    demodu_m(v+1) = Rx_qd_m;
       v=v+2;
end
p= 1;
for z = 1:length(numeric_vector)/2
    if(downsample_rc_in_m(z)<0 && downsample_rc_qd_m(z)<0)
        Rx_in_m_rc = 1;
        Rx_qd_m_rc = 1;
    elseif(downsample_rc_in_m(z)<0 && downsample_rc_qd_m(z)>0)
        Rx_in_m_rc = 1;
        Rx_qd_m_rc = 0;
    elseif(downsample_rc_in_m(z)>0 && downsample_rc_qd_m(z)>0)
        Rx_in_m_rc = 0;
        Rx_qd_m_rc = 0;
    elseif(downsample_rc_in_m(z)>0 && downsample_rc_qd_m(z)<0)
        Rx_in_m_rc = 0;
        Rx_qd_m_rc = 1;
    end
    demodu_m_rc(p)=Rx_in_m_rc;
    demodu_m_rc(p+1) = Rx_qd_m_rc;
       p=p+2;
end
% figure();
% subplot(2,1,1);
% plot(numeric_vector(1:9000));
% title("orig for rect");
% subplot(2,1,2);
% plot(demodu_m(1:9000));
% title("demodu_mlated for rect");
% figure();
% subplot(2,1,1);
% plot(numeric_vector(1:9000));
% title("orig for rc");
% subplot(2,1,2);
% stem(demodu_m_rc(1:9000));
% title("demodu_mlated for rc");
% 
figure();
plot(downsample_rect_qd,downsample_rect_in);
title("constellation of demod rect (sig with memless)");
figure();
plot(downsample_rc_qd,downsample_rc_in);
title("constellation of demod rc(sig with memless)");
figure();
plot(downsample_rect_qd_m,downsample_rect_in_m);
title("constellation of demod rect (sig with memory)");
figure();
plot(downsample_rc_qd_m,downsample_rc_in_m);
title("constellation of demod rc(sig with memory)");
% char_array = num2str(demodu_rc); % Convert numeric array to character array
% char_array = strrep(char_array.', ' ', ''); % Remove any spaces
% binary_matrix = reshape(char_array.', [], 16); % Reshape back to matrix
% audio_integers = bin2dec(binary_matrix); % Convert binary to decimal
% audio_reconstructed = typecast(uint16(audio_integers), 'int16'); % Typecast to int16
% audio_reconstructed_normalized = double(audio_reconstructed) / 32767; % Normalize to [-1, 1]
% % Save the reconstructed audio
% sound(audio_reconstructed_normalized, fs); % Fs is the sampling frequency
char_array_m = num2str(demodu_m_rc); % Convert numeric array to character array
char_array_m = strrep(char_array_m.', ' ', ''); % Remove any spaces
binary_matrix_m = reshape(char_array_m.', [], 16); % Reshape back to matrix
audio_integers_m = bin2dec(binary_matrix_m); % Convert binary to decimal
audio_reconstructed_m = typecast(uint16(audio_integers_m), 'int16'); % Typecast to int16
audio_reconstructed_normalized_m = double(audio_reconstructed_m) / 32767; % Normalize to [-1, 1]
% Save the reconstructed audio
sound(audio_reconstructed_normalized_m, fs); % Fs is the sampling frequency


	
