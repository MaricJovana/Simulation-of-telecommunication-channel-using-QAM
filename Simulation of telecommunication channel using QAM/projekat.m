close all; clear; clc;

%Definisanje frekvencije odabiranja i duzine vremenske ose, tj. broja odbiraka Ns
Fs =384000;
fm = 750;
Ns = (Fs * 128)/fm;
t = (0:1/Fs:(Ns - 1)/Fs);
Nfft=[1024 2048 4096];
Bpou = 3000; 
snr=[5 10 15 20 25];

%frekvencijski spektar
f_2pi=(0:Fs/Nfft(1):Fs-Fs/Nfft(1)); %prikaz od 0 do 2pi
f_pi=(-Fs/2:Fs/Nfft(1):Fs/2-Fs/Nfft(1)); %prikaz od -pi do pi

%generisanje grane u fazi (u odbircima)
Ui_t= cos(2*pi*750*t) + 0.6*cos(2*pi*1500*t) + 0.3*cos(2*pi*2250*t) + 0.1*cos(2*pi*3000*t);
figure(1),subplot(2,1,1), stem(t, Ui_t); title('grana u fazi'); ylim([-5 5]);xlabel('t [µs]');ylabel('ui(t)');


%generisanje grane u kvadraturi (u odbircima)
Uq_t = 4*cos(2*pi*2000*t);
subplot(2,1,2), stem(t, Uq_t); title('grana u kvadraturi');ylim([-5 5]);xlabel('t [µs]');ylabel('uq(t)');

proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs, Ui_t,Uq_t, 'ulazni signali'); 
%mnozenje lokalnim nosiocem, modulacija
f0=Fs/4;
Ui_mod = Ui_t.*cos(2*pi*f0*t);
Uq_mod = Uq_t.*sin(2*pi*f0*t);
Uiq = Ui_mod - Uq_mod; %modulisani signal

proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs, Ui_mod,Uq_mod, 'modulisani signali'); 

%dodavanje aditivnog belog gausovog suma
Uiq_sum = awgn(Uiq, snr(5)-10*log10(Fs/Bpou));

%projektovanje filtra PO
h_po_x1=fir1(50,[f0 f0+3100]/(Fs)); %filtar PO sirine spektra signala B
h_po_x2=fir1(50,[f0-1600 f0+4600]/(Fs)); %filtar PO duplo vece sirine od sirine spektra signala 2B
[H_po_x1,w_po_x1]=freqz(h_po_x1,1);
[H_po_x2,w_po_x2]=freqz(h_po_x2,1);
h_po=[h_po_x1 h_po_x2];

%projektovanje NF filtra
h_nf=fir1(50,3100/Fs);
[H_nf,w_nf]=freqz(h_nf,1);

%filtriranje signala propusnikom opsega
Uiq_po=filter(h_po(1),1, Uiq_sum); 

%idealna demodulacija
Ii_t = Uiq_po*2.* cos(2*pi*f0*t);
Iq_t = Uiq_po*2.* sin(2*pi*f0*t);
proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs, Ii_t,Iq_t, 'idealna demodulacija'); 
Ii_nf=filter(h_nf,1,Ii_t);
Iq_nf=filter(h_nf,1,Iq_t);
proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs,Ii_nf,Iq_nf, 'idealna demodulacija posle NF filtra');

%demodulacija sa faznom greskom; menjati indeks phase_error kako bi se simulacija izvrsila za razlicite vrednosti fazne greske
phase_error = [pi/18 pi/4 pi/2];
Ii_t = Uiq_po*2.* cos(2*pi*f0*t + phase_error(3));
Iq_t = Uiq_po*2.* sin(2*pi*f0*t + phase_error(3));
proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs, Ii_t,Iq_t, 'demodulacija sa faznom greskom'); 
Ii_nf=filter(h_nf,1,Ii_t);
Iq_nf=filter(h_nf,1,Iq_t);
proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs,Ii_nf,Iq_nf,'demodulacija sa faznom greskom posle NF filtra');

%demodulacija sa frekvencijskom greskom ;menjati indeks frequency_error kako bi se simulacija izvrsila za razlicite vrednosti frekvencijske greske
frequency_error = [-0.5 0.5]*1000;
Ii_t = Uiq_po*2.* cos(2*pi*(f0+frequency_error(2))*t);
Iq_t = Uiq_po*2.* sin(2*pi*(f0+frequency_error(2))*t);
proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs, Ii_t,Iq_t,'demodulacija sa frekvencijskom greskom'); 
Ii_nf=filter(h_nf,1,Ii_t);
Iq_nf=filter(h_nf,1,Iq_t);
proracun_i_prikaz_spektra_i_SGSS(Nfft(1),Fs,Ii_nf,Iq_nf,'demodulacija sa frekvencijskom greskom posle NF filtra');