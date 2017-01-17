
%------------------------------------------%
%   NMDOS - projekt                        %
%   Detekcija akorada putem CWT analize    %
%------------------------------------------%
%   Krešimir Rièkoviæ                      %
%   FER, Sijeèanj, 2017.                   %
%------------------------------------------%

clear all
close all

%% Uèitavanje signala
display('Uèitavanje signala');

wavelet = 'cmor32-2';
[signal, fs] = audioread('gitara.wav');
x = (signal(:,1) + signal(:,2));
sub = 32;
 
x = x(1:sub:length(x),1);   % poduzorkovanje faktorom 32
T = 1/fs*sub;               % period
t = [0:length(x)-T]' * T;   % vektor vremena

%% Wavelet analiza
display('Wavelet analiza');

[psi,xval] = wavefun(wavelet);
 
a3 = scal2frq(1, wavelet, T) / 329.6;            % E4 ton (329.6 Hz) - skala

figure(1)
subplot(2,1,1),plot(xval,real(psi));
title('Realni dio');
subplot(2,1,2),plot(xval,imag(psi));
title('Imaginarni dio');
% saveas(gcf,'Valiæ.png');
  
%              F#,  F,  E, D#,  D, C#,  C,  B, A#,  A, G#, G
%        2.^([ -9, -8, -7, -6, -5, -4, -3, -2, -1, 0 , 1 , 2 ]/12) * frek_a

step = 0.0625;
octave = a3 * 2.^([-9:step:2+(1/step-1)*step]/12);
 
scales = [octave 2*octave 4*octave];            % skale
frequencies = scal2frq(scales, wavelet, T);     % frekvencije
 
% CWT
X = cwt(x, scales, wavelet);

% vizualizacija
figure(2), h=imagesc(t, frequencies, abs(X)); ylabel('Hz'), xlabel('s'),
set(get(h, 'Parent'), 'YDir', 'normal')
set(get(h, 'Parent'), 'YScale', 'log')

oct = a3 * 2.^([-9:2]/12);
fre = fliplr(scal2frq([oct 2*oct 4*oct], wavelet, T));
set(get(h, 'Parent'), 'YTick', fre)
grid
grid minor
% saveas(gcf,'T-F graf.png');

%% Analiziranje akordiju
display('Analiziranje akordiju');

lower_thresh = 0.1;
upper_thresh = 0.2;
position = zeros(10,length(X(1,:)));       % 5 pozicija (poèetak, kraj) : x-os
energy = zeros(10,length(X(1,:)));         % 5 pozicija (energija, pozicija) : x-os
tones(length(X(1,:)),6) = char(0);

fileID = fopen('tonovi.txt','w');

% dur stringovi

C  = 'C   CEGCECGEC';
d  = 'C#  dFadFdaFd';
D  = 'D   DgADgDAgD';
e  = 'D#  eGbeGebGe';
E  = 'E   EaBEaEBaE';
F  = 'F   FACFAFCAF';
g  = 'F#  gbdgbgdbg';
G  = 'G   GBDGBGDBG';
a  = 'G#  aCeaCaeCa';
A  = 'A   AdEAdAEdA';
b  = 'A#  bDFbDbFDb';
B  = 'B   BegBeBgeB';

% mol stringovi

Cm = 'Cm  CeGCeCGeC';
dm = 'C#m dEadEdaEd';
Dm = 'Dm  DFADFDAFD';
em = 'D#m egbegebge';
Em = 'Em  EGBEGEBGE';
Fm = 'Fm  FaCFaFCaF';
gm = 'F#m gAdgAgdAg';
Gm = 'Gm  GbDGbGDbG';
am = 'G#m aBeaBaeBa';
Am = 'Am  ACEACAECA';
bm = 'A#m bdFbdbFdb';
Bm = 'Bm  BDgBDBgDB';

chords = [C; d; D; e; E; F; g; G; a; A; b; B; Cm; dm; Dm; em; Em; Fm; gm; Gm; am; Am; bm; Bm];

for i=1:length(X(1,:))
    
    joint_flag = false;     % zastavica za karakteriziranje podruèja kao jedan ton
    tone_count = 0;         % brojaè tonova
    max_E = 0;              % maksimalna energija - daje informaciju o frekvenciji tona
    joint_count = 0;        % brojaè za odbacivanje šumova koji su debljine do 3
    
    % Segment: 0.000070 seconds (bez tona), 0.001040 seconds (sa tonom)
    for j=1:length(X(:,1))
        
        if(joint_flag == false)
            if(abs(X(j,i)) > upper_thresh)
                joint_flag = true;
                position(tone_count*2+1,i) = j;               % oznaèi poèetak tranzicije u ton
                joint_count = joint_count + 1;
            end
            
        else
            if(abs(X(j,i)) < lower_thresh)
                joint_flag = false;
                max_E = 0;
                
                if(joint_count > 3)
                    position(tone_count*2+2,i) = j;               % oznaèi kraj tranzicije u ton
                    tone_count = tone_count + 1;
                else
                    position(tone_count*2+1,i) = 0;
                    energy(tone_count*2+1,i) = 0;
                    energy(tone_count*2+2,i) = 0;
                end 
                joint_count = 0;
                
            elseif(abs(X(j,i)) > max_E)
                max_E = abs(X(j,i));
                energy(tone_count*2+1,i) = max_E;
                energy(tone_count*2+2,i) = frequencies(j);
                joint_count = joint_count + 1;
            end
        end
    end
    
    % Segment: 0.000000 seconds (bez tona), 0.000797 seconds (s tonom)
    for j=1:tone_count
        tone = energy(j*2,i);
        for k=0:2
            if(tone >= 127.15*2^k && tone < 134.7*2^k)
                tones(i,j) = strcat('C');
            elseif(tone >= 134.7*2^k && tone < 142.71*2^k)
                tones(i,j) = strcat('d');
            elseif(tone >= 142.17*2^k && tone < 151.2*2^k)
                tones(i,j) = strcat('D');
            elseif(tone >= 151.2*2^k && tone < 160.19*2^k)
                tones(i,j) = strcat('e');
            elseif(tone >= 80.1*2^k && tone < 84.86*2^k)
                tones(i,j) = strcat('E');
            elseif(tone >= 84.86*2^k && tone < 89.91*2^k)
                tones(i,j) = strcat('F');
            elseif(tone >= 89.91*2^k && tone < 95.25*2^k)
                tones(i,j) = strcat('g');
            elseif(tone >= 95.25*2^k && tone < 100.92*2^k)
                tones(i,j) = strcat('G');
            elseif(tone >= 100.92*2^k && tone < 106.92*2^k)
                tones(i,j) = strcat('a');
            elseif(tone >= 106.92*2^k && tone < 113.27*2^k)
                tones(i,j) = strcat('A');
            elseif(tone >= 113.27*2^k && tone < 120.01*2^k)
                tones(i,j) = strcat('b');
            elseif(tone >= 120.01*2^k && tone < 127.15*2^k)
                tones(i,j) = strcat('B');
            end
        end
    end
    
    % Segment: 0.000266 seconds
    for j=1:length(tones(i,:))
        for k=j+1:length(tones(i,:))
            if(tones(i,j) == tones(i,k))
                tones(i,k:end) = [tones(i,k+1:end) blanks(1)];
                k = k - 1;
            end
        end
    end
    
    if(~istone(tones(i,3)))
        tones(i,:) = blanks(length(tones(i,:)));
    end
    
    k = 0;       % brojaè tonova u stringu
    while(istone(tones(i,k+1)))
        k = k + 1;
    end

    chord_exists = false;
    for j=1:length(chords)
        if(~isempty(strfind(chords(j,:),tones(i,1:k))))
            fprintf(fileID,chords(j,1:k));
            fprintf(fileID,'\n');
            chord_exists = true;
        end
    end
    if(~chord_exists)
        fprintf(fileID,'\n');
    end
end

fclose(fileID);
