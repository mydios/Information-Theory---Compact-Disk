% the purpose of this script is to check whether dec(enc(x)) == x
% same for the delay and interleave functions.

% make cd
Fs = 22.05e3;
cd = AudioCD(Fs,1,8);

% --------- CIRC_dec_delay_interleave(CIRC_enc_delay_interleave(input)) == input
n_frames = 2;
input = zeros(n_frames*24,1,'uint8');

% make some dummy data
for i = 1:n_frames*24
   input(i) = i; 
end

[output, n_frames_output] = cd.CIRC_enc_delay_interleave(input,n_frames);
erasures = zeros(n_frames_output*24,1,'logical');
[output2, flags_out, n_frames_output2] = cd.CIRC_dec_deinterleave_delay(output, erasures, n_frames_output);

fprintf('CIRC_dec_delay_interleave(CIRC_enc_delay_interleave(input)) == input')
isequal(input, output2)

% --------- CIRC_dec_C2(CIRC_enc_C2(input)) == input
n_frames = 2;
input = zeros(n_frames*24,1,'uint8');

% make some dummy data
for i = 1:n_frames*24
   input(i) = i; 
end

[output, n_frames_output] = cd.CIRC_enc_C2(input,n_frames);
erasures = zeros(n_frames_output*28,1,'logical');
[output2, flags_out, n_frames_output2] = cd.CIRC_dec_C2(output, erasures, n_frames_output);

fprintf('CIRC_dec_C2(CIRC_enc_C2(input)) == input')
isequal(input, output2)

% -------- CIRC_dec_delay_unequal(CIRC_enc_delay_unequal(input)) == input
n_frames = 2;
input = zeros(n_frames*28,1,'uint8');

% make some dummy data
for i = 1:n_frames*28
   input(i) = i; 
end

[output, n_frames_output] = cd.CIRC_enc_delay_unequal(input,n_frames);
erasures = zeros(n_frames_output*28,1,'logical');
[output2, flags_out, n_frames_output2] = cd.CIRC_dec_delay_unequal(output, erasures, n_frames_output);

fprintf('CIRC_dec_delay_unequal(CIRC_enc_delay_unequal(input)) == input')
isequal(input, output2)

% -------- CIRC_dec_C1(CIRC_enc_C1(input)) == input
n_frames = 2;
input = zeros(n_frames*28,1,'uint8');

% make some dummy data
for i = 1:n_frames*28
   input(i) = i; 
end

[output, n_frames_output] = cd.CIRC_enc_C1(input,n_frames);
[output2, flags_out, n_frames_output2] = cd.CIRC_dec_C1(output, n_frames_output);

fprintf('CIRC_dec_C1(CIRC_enc_C1(input)) == input')
isequal(input, output2)

% ------- CIRC_dec_delay_inv(CIRC_enc_delay_inv(input)) == input
n_frames = 2;
input = zeros(n_frames*32,1,'uint8');

% make some dummy data
for i = 1:n_frames*32
   input(i) = i; 
end

[output, n_frames_output] = cd.CIRC_enc_delay_inv(input,n_frames);
[output2, n_frames_output2] = cd.CIRC_dec_delay_inv(output, n_frames_output);

fprintf('CIRC_dec_delay_inv(CIRC_enc_delay_inv(input)) == input')
isequal(input, output2)
