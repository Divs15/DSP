%LAB4
%%LAB 3.1.1
%% load data 
% load DNA sequence we want to work with
hbb = genbankread('hbb_region_chr11(1).gb');
% display the method we are using
disp('Non-Parallelized method:')
tic
% run STFT function (provided on BBLEARN)
Threebaseperiodicity_vs_position_old=threebasefreq_stft(hbb.Sequence,1000,1024);
%%
%{
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

%}
% use toc to mark an end time
% tic & toc will print the elapsed time for commands between tic and toc
% in our case, it print out how long it took to run STFT function
toc
% plot the result
figure(1)
plot(Threebaseperiodicity_vs_position_old)
title('Threebaseperiodicity_vs_position_old')

Threebaseperiodicity_vs_position = {};
worker_num = 2;
% set parpool (check the Matlab help)
delete(gcp('nocreate'))
parpool(worker_num);
disp('Parallelized method:')
tic
parfor i=1:2
    Threebaseperiodicity_vs_position{i} = p_threebasefreq_stft(hbb.Sequence, 1000, 1024, 4, i);
end
toc
%% Function p_threebasefreq_stft 
%{
function Threebaseperiodicity_vs_position = p_threebasefreq_stft(Sequence, WINDOW_LENGTH, NFFT, division, dindex)
whole_length = length(Sequence);
seq_div = [];
step_size = floor(whole_length/division);
start_point = (dindex-1)*step_size+1;
if dindex == division
    end_point = whole_length;
else
    end_point = dindex*step_size+WINDOW_LENGTH-1;
end
seq_div = Sequence(start_point:end_point);
Threebaseperiodicity_vs_position = threebasefreq_stft (seq_div, WINDOW_LENGTH, NFFT);
end
    
%}
%%
%LAB3.2.1
Threebaseperiodicity_vs_position_New = [];
for i=1:2
    
    Threebaseperiodicity_vs_position_New = [Threebaseperiodicity_vs_position_New,Threebaseperiodicity_vs_position(i)];
end
A = cell2mat(Threebaseperiodicity_vs_position_New);
figure(2)
plot(A)
title('Threebaseperiodicity_vs_position_New')
%% check the difference between old and new values:
diff = sum(abs(A - Threebaseperiodicity_vs_position_old));
disp('Absolute Sum Error:')
disp(diff)
