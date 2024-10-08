%Author:    Dimme de Groot
%Date:      July 2024
%Descr:     This code illustrates how to use the microphone spotformer object when precomputing the weights.

clear all
close all

listen_flag = true;     %Set to false if you do not want to listen to the audio!

addpath RIR-generator/  %Path to RIR generator by E. Habets (see license)
addpath clenquad/       %Path to Clenshaw-Curtis quadrature by G. von Winckel (see license)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the received audio (audioRec). This is the mixture of the target audio (audioRecTAR) and the interfering (audioRecINT) audio.  %
% AudioPlayClean is the input signal. NN is the microphone nearest to the user                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) First set the microphone, interferer and target locations

% (1a) Microphone (receiver) positions. Here they are placed in a circular array of 10 cm radius. 
Nmic = 8;                               % Number of microphones
theta = linspace(0, 2*pi, Nmic+1);      %[rad], simulate circular mic array 
theta = theta(1:Nmic).';                            
r = 0.1;                                %[m], radius circular array
loc_mic = [r*cos(theta), r*sin(theta), zeros(Nmic,1)]; 
loc_mic = loc_mic + [2, 2, 1];          %[m], microphone positions

% (1b) Loudspeaker (interferer) position
loc_loud = [6.010, 2.019, 1.175];       %[m], FL loudspeaker position
            
% (1c) Person (target) position
loc_pers = [6,4, 1.100];        %[m], person position (target)

% (1d) Plot the problem setup!
plotLayout(loc_loud, loc_pers, loc_mic)
drawnow()

% (2) Compute the audio received at the microphones 
noiseStrength = 0.3;    %Increase/decrease to add/remove noise
sound_vel = 342;        %[m/s], speed of sound
fs = 16000;             %[Hz], sample frequency
[audioRec, audioRecTAR, audioRecINT, audioPlayClean, NN] = fnc_computeReceivedAudio(loc_loud, loc_mic, loc_pers, noiseStrength, sound_vel, fs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We are now gonna set the spotformer object, compute the microphone weights and compute the output signal % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (1) Inputs related to integration
N_int = 10;             %[-], number of integration points
IntWinRad = 0.2;        %[m], interferer window radius (3 sigma)              
TarWinRad = 0.5;        %[m], target window radius (3 sigma)              

% (2) Other inputs: processing
flag_full_axis = false; %[-], True for full frequency axis [-Fs/2, Fs/2). False for [0, Fs/2]. True is not tested

t_frame = 0.016;        %[s], analysis window length
t_pad = 0.016;          %[s], padding window length

rebRatio = 0.0;         %Term describing direct to reverberant component (I set this number arbitrarily)
numSigma2 = 10^-9;      %Term for dealing with numerical inaccuracies stemming from e.g. numerical integration. Effectively regularises the result by enforcing positive definiteness
nSigma2 = 0;            %Term which can be set in case of microphone self noise.

% (3) Analysis and synthesis window
analysis_window = "sqrthann";
synthesis_window = "sqrthann";

% (4) Create microphone spotformer object
MicSpot = MicSpotformer(sound_vel, fs,t_frame, t_pad, N_int, IntWinRad, TarWinRad, nSigma2, numSigma2, rebRatio, flag_full_axis, analysis_window, synthesis_window);

% (5) Compute spotformer weights
MicSpot.fnc_comp_weights(loc_loud, loc_pers, loc_mic)

% (6) Compute output spotformer
output = MicSpot.comp_output(audioRec);

%%%%%%%%%%%%%%%%%%%
% Listen to audio %
%%%%%%%%%%%%%%%%%%%
if listen_flag == true
    disp('playing audio...')
    disp('Clean audio at nearest microphone...')
    soundsc(audioRecTAR(:, NN), fs)
    pause(length(audioRecTAR(:,NN))/fs)

    disp('Noisy mixture at nearest microphone...')
    soundsc(audioRec(:, NN), fs)
    pause(length(audioRecTAR(:,NN))/fs)

    disp('Output microphone spotformer...')
    soundsc(output, fs)
end

rmpath clenquad/
rmpath RIR-generator/

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some extra functions for generating the received audio and plotting the microphones and loudspeakers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotLayout(locInterferer, locTarget, locReceiver)
    figure
    hold on
    scatter3(locInterferer(:,1), locInterferer(:,2), locInterferer(:,3), 'filled');
    scatter3(locTarget(:,1), locTarget(:,2), locTarget(:,3), 'filled');
    scatter3(locReceiver(:,1), locReceiver(:,2), locReceiver(:,3), 'filled');
    axis tight
    axis equal
    grid on
    legend('Interferer', 'Target', 'Microphones')
end

function [audioRec, audioRecTAR, audioRecINT, audioTar, NN] = fnc_computeReceivedAudio(locInterferer, locReceiver, locTarget, sigma2_noise, c, Fs)
    if nargin == 3
        c = 342;        %[m/s], sound velocity
        Fs = 16000;     %[Hz],  sampling frequency
        sigma2_noise = 0.3; %parameter setting the noise intensity
    end

    if min(locInterferer, [], 'all') < 0.2 || min(locReceiver, [], 'all') < 0.2 || min(locTarget, [], 'all') < 0.2
        disp("Please ensure all coordinates are larger than 0.2 m. At least one of your coordinates does not adhere to this.")
    end

    %parametrize RIR
    L = [7 6 3];                %[m], Room dimensions [Lx, Ly, Lz] 
    beta = 0;                   %[s], reverberation time
    nRir = 4096;                %[-], Number of samples in RIR

    if beta ~= 0
        aS_x = L(2)*L(3)*(1-beta(1)^2+1-beta(2)^2);
        aS_y = L(1)*L(3)*(1-beta(3)^2+1-beta(4)^2);
        aS_z = L(1)*L(2)*(1-beta(5)^2+1-beta(6)^2);
        V = L(1)*L(2)*L(3);
        T60 = 24*log(10)*V/(c*(aS_x+aS_y+aS_z))
    end

    %Some handy numbers 
    N_int = size(locInterferer,1);
    N_tar = size(locTarget, 1); 
    N_rec = size(locReceiver,1);

    if N_tar == 1
        audioTar = audioread('Excerpt/sample_16kHz.wav');
    elseif N_tar == 2
        audioTar(:,1) = audioread('Excerpt/sample_16kHz.wav');
        audioTar(:,2) = zeros(size(audioTar(:,1)));
        train = load('train');
        audioTarTMP = resample(train.y, Fs, train.Fs);
        audioTar(1:length(audioTarTMP),2) = audioTarTMP;
    else
        disp("Some modifications of fnc_computeReceivedAudio are required for more than two persons :S");
    end
    audioTar = [zeros(5*Fs,N_tar); audioTar; zeros(3*Fs,N_tar)];
    audioInt = sigma2_noise*randn(length(audioTar), N_int);

    for i=1:N_rec
        for j=1:N_int
            RIR = rir_generator(c, Fs, locReceiver(i,:), locInterferer(j,:), L, beta, nRir);
            try
                audioRecINT(:,i) = audioRecINT(:,i) + conv(RIR,audioInt(:,j));
            catch 
                audioRecINT(:,i) = conv(RIR,audioInt(:,j));
            end
        end
        for j=1:N_tar
            RIR = rir_generator(c, Fs, locReceiver(i,:), locTarget(j,:), L, beta, nRir);
            try
                audioRecTAR(:,i) = audioRecTAR(:,i) + conv(RIR,audioTar(:,j));
            catch 
                audioRecTAR(:,i) = conv(RIR,audioTar(:,j));
            end
        end
    end
    audioRec = audioRecINT+audioRecTAR;

    %Give estimate SNR
    [~, NN] = min(vecnorm(locReceiver-locTarget(1,:),2,2));
    S = audioRecTAR(:,NN);
    N = audioRecINT(:,NN);
    INDX = find(abs(S)>0.01*max(abs(S)));
    S = S(INDX);
    N = N(INDX);
    SNR = 20*log10(norm(S)/norm(N));
    disp("estimated SNR at reference microphone (" + num2str(NN) + ") is " + num2str(SNR) + " dB") 
end