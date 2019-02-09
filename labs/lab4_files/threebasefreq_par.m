function output = threebasefreq_par(DNA_SEQUENCE, WINDOW_LENGTH, NFFT, division, dindex)
break_points = (length(DNA_SEQUENCE)-WINDOW_LENGTH)/division;

seq = DNA_SEQUENCE((dindex+break_points*(dindex-1)):break_points*dindex+WINDOW_LENGTH);

DNA_len = length(seq);
output=zeros(1,DNA_len-WINDOW_LENGTH+1);
coding = seq; 
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
    output(i)=STFT(floor(NFFT/3));
end