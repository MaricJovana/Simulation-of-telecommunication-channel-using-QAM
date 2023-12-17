function proracun_i_prikaz_spektra_i_SGSS(Nfft,Fs,Ui_t,Uq_t, recenica)

f_fja=(0:Fs/Nfft:Fs-Fs/Nfft); %f osa
figure,
subplot(2,1,1), stem( f_fja, abs(fft(Ui_t,Nfft))); title('grana u fazi ', recenica); xlim([0 Fs]);xlabel('f [Hz]');
subplot(2,1,2), stem( f_fja, abs(fft(Uq_t,Nfft))); title('grana u kvadraturi', recenica); xlim([0 Fs]);xlabel('f [Hz]');
sgtitle('frekvencijski spektar');

Ui_f = zeros((65536/Nfft), Nfft);
Uq_f = zeros((65536/Nfft), Nfft);
Ui_f(1, :) = abs(fft(Ui_t(1:Nfft),Nfft));
Uq_f(1, :) = abs(fft(Uq_t(1:Nfft),Nfft));
for br =2:(65536/Nfft)
    Ui_f(br,:) = abs(fft(Ui_t((br-1)*Nfft:br*Nfft),Nfft));
    Uq_f(br,:) = abs(fft(Uq_t((br-1)*Nfft:br*Nfft),Nfft));
end
%proracun SGSS signala
Ui_zbir_2 = sum((Ui_f.^2)./(Fs*Nfft));
Uq_zbir_2 = sum((Uq_f.^2)./(Fs*Nfft));
%SGSS po formuli sa postavke zadatka
Ui_usrednjeno_2 = Ui_zbir_2./(65536/Nfft); 
Uq_usrednjeno_2 = Uq_zbir_2./(65536/Nfft);

%prikaz SGSS
figure,
subplot(2,1,1), stem( f_fja, Ui_usrednjeno_2); title('grana u fazi usrednjena', recenica); xlim([0 Fs]);xlabel('f [Hz]');ylabel('sgss');
subplot(2,1,2), stem( f_fja, Uq_usrednjeno_2); title('grana u kvadraturi usrednjena', recenica); xlim([0 Fs]);xlabel('f [Hz]');ylabel('sgss');
sgtitle('SGSS');

end