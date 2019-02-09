function Threebaseperiodicity_vs_position = threebasefreq_stft (DNA_SEQUENCE, WINDOW_LENGTH, NFFT)
DNA_len = length(DNA_SEQUENCE);
Threebaseperiodicity_vs_position=zeros(1,DNA_len-WINDOW_LENGTH+1);
coding = DNA_SEQUENCE; 
coding_A = (upper(coding)=='A'); % find A bases and set them to 1
coding_T = (upper(coding)=='T'); % find T bases and set them to 1
coding_G = (upper(coding)=='G'); % find G bases and set them to 1
coding_C = (upper(coding)=='C'); % find C bases and set them to 1
for i = 1:DNA_len-WINDOW_LENGTH+1
    stc_A = coding_A(i:i+WINDOW_LENGTH-1);
    stc_T = coding_T(i:i+WINDOW_LENGTH-1);
    stc_G = coding_G(i:i+WINDOW_LENGTH-1);
    stc_C = coding_C(i:i+WINDOW_LENGTH-1);
    STFT = abs(fft(stc_A,NFFT)).^2+abs(fft(stc_T,NFFT)).^2 ...
    +abs(fft(stc_G,NFFT)).^2+abs(fft(stc_C,NFFT)).^2; % FFT of the sequence
    Threebaseperiodicity_vs_position(i)=STFT(floor(NFFT/3));
end