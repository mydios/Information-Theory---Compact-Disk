classdef AudioCD
    % AudioCD class: simulation of the CIRC encoding/decoding and interpolation of the CD digital audio standard
    %
    % Author: Johannes Van Wonterghem, Jan 2017
    %
    % See the static test() method for an example usage of this class
    
    properties
        Fs; % Sample rate of the audio
        
        configuration; % % 0: no CIRC; 1: CIRC as described in standard; 2: Concatenated RS, no interleaving; 3: Single 32,24 RS
        max_interpolation; % The maximum number of interpolated audio samples
        
        enc2; % Outer RS code of the CIRC
        dec2;
        gpoly2;
        
        enc1; % Inner RS code of the CIRC
        dec1;
        gpoly1;
        
        gpoly_8_parity; % Single RS code for configuration 3
        enc_8_parity;
        dec_8_parity;
        
        cd_bits; % Bits written to disk (before EFM)
        
        scaled_quantized_padded_original; % Reference to compare the output of readCD to
    end
    
    methods
        
        function obj=AudioCD(Fs,configuration,max_interpolation)
            % Constructor of the AudioCD class
            % INPUT:
            % -Fs: The sample rate of the audio
            % -configuration: 0: no CIRC; 1: CIRC as described in standard; 2: Concatenated RS, no interleaving; 3: Single 32,24 RS
            % -max_interpolation: The maximum number of interpolated audio samples
            % OUTPUT:
            % -obj: the AudioCD object
            
            obj.Fs = Fs;
            obj.max_interpolation = max_interpolation;
            
            % Initialize the RS encoders and decoders
            primpoly = [1,0,0,0,1,1,1,0,1]; % Primitive polynomial from the standard
            if (configuration == 1) || (configuration == 2)
                obj.gpoly2 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc2 = comm.RSEncoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec2 = comm.RSDecoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly,'ErasuresInputPort',true);
                
                obj.gpoly1 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc1 = comm.RSEncoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec1 = comm.RSDecoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly,'ErasuresInputPort',true);
            elseif configuration == 3
                obj.gpoly_8_parity = rsgenpoly(255,247,bi2de(fliplr(primpoly)),0);
                obj.enc_8_parity = comm.RSEncoder(255,247,obj.gpoly_8_parity,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec_8_parity = comm.RSDecoder(255,247,obj.gpoly_8_parity,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
            end
            obj.configuration = configuration;
        end
        
        function obj=writeCd(obj,audiofile)
            % Write an audiofile to the CD
            % INPUT:
            % -obj: the AudioCD object
            % -audiofile: [Nsamples x 2] matrix containing the left and right audio track as samples of the double datatype
            % OUTPUT:
            % -obj: the updated AudioCD object
            
            assert(size(audiofile,2) == 2);
            
            xscaled = audiofile / max(max(abs(audiofile))); % normalize to -1:1
            x = uencode(xscaled,16); % convert to 16 bit signed values
            
            xlr16 = reshape(x',[],1); % serialize left and right audio channel
            xlr8 = typecast(xlr16,'uint8'); % split into 8 bit words
            
            xlr8_padded = [xlr8 ; zeros(24-(rem(numel(xlr8)-1,24)+1),1)]; % pad with zeros to fill an integer number of frames
            n_frames = numel(xlr8_padded)/24; % every frame contains 24 8 bit words
            
            ylr16 = typecast(uint8(xlr8_padded),'uint16');
            y = reshape(ylr16,2,[])';
            obj.scaled_quantized_padded_original = udecode(y,16); % Reference to compare the output of readCD to
            
            switch obj.configuration
                case 0 % no CIRC
                    xlrb = de2bi(xlr8,8);
                case 1 % CIRC as described in standard
                    [delay_interleaved,n_frames] = obj.CIRC_enc_delay_interleave(xlr8_padded,n_frames);
                    [C2_encoded,n_frames] = obj.CIRC_enc_C2(delay_interleaved,n_frames);
                    [delay_unequal,n_frames] = obj.CIRC_enc_delay_unequal(C2_encoded,n_frames);
                    [C1_encoded,n_frames] = obj.CIRC_enc_C1(delay_unequal,n_frames);
                    [delay_inv,n_frames] = obj.CIRC_enc_delay_inv(C1_encoded,n_frames);
                    
                    xlrb = de2bi(delay_inv,8);
                case 2 % Concatenated RS, no interleaving
                    [C2_encoded,n_frames] = obj.CIRC_enc_C2(xlr8_padded,n_frames);
                    [C1_encoded,n_frames] = obj.CIRC_enc_C1(C2_encoded,n_frames);
                    
                    xlrb = de2bi(C1_encoded,8);
                case 3 % Single 32,24 RS
                    [encoded,n_frames] = obj.C3_enc_8_parity(xlr8_padded,n_frames);
                    
                    xlrb = de2bi(encoded,8);
                otherwise
                    error('Invalid configuration selected');
            end
            
            xlrbserial = reshape(xlrb',[],1);
            
            obj.cd_bits = logical(xlrbserial);
        end
        
        function obj=bitErrorsCd(obj,p)
            % Add uniform bit errors to cd
            % INPUT:
            % -obj: the AudioCD object
            % -p: the bit error probability, i.e., an obj.cd_bits bit is flipped with probability p
            % OUTPUT:
            % -obj: the updated AudioCD object
            
            noise = rand(size(obj.cd_bits))<p;
            
            obj.cd_bits = xor(obj.cd_bits,noise);
        end
        
        function obj=scratchCd(obj,length_scratch,location_scratch)
            % Add a scratch to the cd
            % INPUT:
            % -obj: the AudioCD object
            % -length_scratch: the length of the scratch (in number of bit)
            % -location_scratch: the location of the scratch (in bit offset from start of obj.cd_bits)
            % OUTPUT:
            % -obj: the updated AudioCD object
            
            obj.cd_bits(location_scratch:min(location_scratch+length_scratch-1,numel(obj.cd_bits))) = 0;
        end
        
        function [audio_out,interpolation_flags]=readCd(obj)
            % Read an audiofile from the CD
            % INPUT:
            % -obj: the AudioCD object, containing the data in obj.cd_bits
            % OUTPUT:
            % -audio_out: [Nsamples x 2] matrix containing the left and right audio track as samples of the double datatype
            % -interpolation_flags: [Nsamples x 2] matrix containing a 0 where no erasure was flagged, a 1 where an erasure was interpolated and a -1
            % where interpolation failed
            
            
            ylrb = reshape(obj.cd_bits,8,[])';
            ylr8 = bi2de(ylrb);
            
            switch obj.configuration
                case 0 % no CIRC
                    ylr16 = typecast(uint8(ylr8),'uint16');
                    y = reshape(ylr16,2,[])';
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                case 1
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [delay_inv,n_frames] = obj.CIRC_dec_delay_inv(ylr8,n_frames);
                    [C1_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C1(delay_inv,n_frames);
                    [delay_unequal,erasure_flags,n_frames] = obj.CIRC_dec_delay_unequal(C1_decoded,erasure_flags,n_frames);
                    [C2_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C2(delay_unequal,erasure_flags,n_frames);
                    [deinterleave_delay,erasure_flags,n_frames] = obj.CIRC_dec_deinterleave_delay(C2_decoded,erasure_flags,n_frames);
                    
                    ylr16 = typecast(uint8(deinterleave_delay),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                    
                case 2 % Concatenated RS, no interleaving
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [C1_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C1(ylr8,n_frames);
                    erasure_flags_t = erasure_flags;
                    [C2_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C2(C1_decoded,erasure_flags,n_frames);
                    
                    if(numel(erasure_flags)  ~= numel(C2_decoded))
                        disp('Something wrong!');
                    end
                    
                    ylr16 = typecast(uint8(C2_decoded),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                case 3 % Single 32,24 RS
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [decoded,erasure_flags,n_frames] = obj.C3_dec_8_parity(ylr8,n_frames);
                    
                    ylr16 = typecast(uint8(decoded),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                otherwise
                    error('Invalid configuration selected');
            end
            
        end
        
        function [output,n_frames] = CIRC_enc_delay_interleave(~,input,n_frames)
            % CIRC Encoder: Delay of 2 frames + interleaving sequence
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % There are 24 symbols in each frame which make up 12 words.
            % Each word can at max be delayed up to 2 frames.
            % As there are n_frames in total, the words in the last frame
            % can be put in the n_frames+2'th frame
            % Initialize an array that can fit symbols of all the frames
            output = zeros((n_frames+2)*24,1,'uint8');
            
            % add delay + do interleaving
            for i=0:n_frames-1
                
                frame_index = (i)*24;
                delayed_frame_index = (i+2)*24;
                
                % Retrieve the corresponding frame from the input
                frame = input(frame_index + 1:frame_index + 24);
                
                % All of the even words get delayed by 2 frames and put at
                % the start of the frame
                output(delayed_frame_index + 1 :delayed_frame_index+2) = frame(1:2); %L6n+0
                output(delayed_frame_index + 3:delayed_frame_index+4) = frame(9:10); %L6n+2
                output(delayed_frame_index + 5:delayed_frame_index+6) = frame(17:18); %L6n+4
                
                output(delayed_frame_index + 7:delayed_frame_index+8) = frame(3:4); %R6n+0
                output(delayed_frame_index + 9:delayed_frame_index+10) = frame(11:12); %R6n+2
                output(delayed_frame_index + 11:delayed_frame_index+12) = frame(19:20); %R6n+4
                
                
                % The uneven words are not delayed but the left and right
                % channel also do get split and put at the bottom of the
                % frame
                output(frame_index + 13:frame_index+14) = frame(5:6); %L6n+1
                output(frame_index + 15:frame_index+16) = frame(13:14); %L6n+3
                output(frame_index + 17:frame_index+18) = frame(21:22); %L6n+5
                
                output(frame_index + 19:frame_index+20) = frame(7:8); %R6n+1
                output(frame_index + 21:frame_index+22) = frame(15:16); %R6n+3
                output(frame_index + 23:frame_index+24) = frame(23:24); %R6n+5
                
            end
            
            % Increase the total number of frames to account for the new
            % frames
            n_frames = n_frames + 2;
            
        end
        
        function [output,n_frames] = CIRC_enc_C2(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C2)
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder
            % -n_frames: the length of the output expressed in frames
            
            % Prepare the output of this block for 28 symbols. 4 parity
            % symbols get added in this step.
            output = zeros(n_frames*28,1,'uint8');
            
            for i = 0:n_frames-1
                
                frame_index = i*24;
                new_frame_index = i*28;
                
                % Select the frame
                frame = input(frame_index+1:frame_index+24);
                
                % Encode the selected frame
                enc_out = step(obj.enc2,frame);
                
                % Place the 4 parity symbols in the middle of the frame
                % The encoder, by default, places these at the end of the
                % frame.
                output(new_frame_index+1:new_frame_index+12) = enc_out(1:12);
                output(new_frame_index+13:new_frame_index+16) = enc_out(25:28);
                output(new_frame_index+17:new_frame_index+28) = enc_out(13:24);
            end
        end
        
        function [output,n_frames] = CIRC_enc_delay_unequal(~,input,n_frames)
            % CIRC Encoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % Prepare output buffer with additional frames
            output = zeros((n_frames+27*4)*28,1,'uint8');
            
            % The first symbol does not get delayed so only 27 symbols need
            % to be adjusted
            for i=0:n_frames-1
                
                frame_index = i*28;
                
                for j=0:27
                    
                    delay = j*4*28;
                    symbol_index = j+1;
                    
                    % Map the input frame on the corresponding delayed output
                    % frame
                    output(frame_index + delay + symbol_index) = input(frame_index + symbol_index);
                end
            end
            n_frames = n_frames + 27*4;
        end
        
        function [output,n_frames] = CIRC_enc_C1(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C1)
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder
            % -n_frames: the length of the output expressed in frames
            
            % Input are 28 symbols / frame, output is 32 symbols / frame
            output = zeros(n_frames*32,1,'uint8');
            
            for i = 0:n_frames-1
                
                frame_index = i*28;
                new_frame_index = i*32;
                
                output(new_frame_index + 1 : new_frame_index + 32) = step(obj.enc1,input(frame_index + 1 : frame_index + 28));
                % No need to correct here like in C2 as the parity bits
                % already get put at the bottom of the frame
            end
        end
        
        function [output,n_frames] = CIRC_enc_delay_inv(~,input,n_frames)
            % CIRC Encoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % The output will grow by one frame as the last symbol will be
            % delayed with 1 frame
            output = zeros((n_frames+1)*32,1,'uint8');
            
            for i=0:n_frames-1
                
                frame_index = i*32;
                delayed_frame_index = (i+1)*32;
                
                for j=0:31
                    
                    symbol_index = j + 1;
                    index = frame_index; % frame_index by default
                    
                    % Use delayed index on even frames
                    if mod(j,2) == 0
                        index = delayed_frame_index;
                    end
  
                    output(index+symbol_index) = input(frame_index+symbol_index);
                        
                    % The middle 4 symbols up until 28 and last 4
                    % symbols need to be inverted
                    if (j >= 12 && j <= 15) || j>=28
                        output(index+symbol_index) = bitcmp(output(index+symbol_index),'uint8');
                    end    

                end
            end
            
            n_frames = n_frames + 1;
            
        end
        
        function [output,n_frames] = CIRC_dec_delay_inv(~,input,n_frames)
            % CIRC Decoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_inv(obj.CIRC_enc_delay_inv(x)) == x!
            
            % Perform the reverse of the enc
            % Just like in the enc version, a 1 frame delay is performed
            % but now on the uneven symbols. This again adds a extra frame
            output = zeros((n_frames+1)*32,1,'uint8');
            
            for i=0:n_frames-1
                
                frame_index = i*32;
                delayed_frame_index = (i+1)*32;
                
                for j=0:31
                    
                    symbol_index = j + 1;
                    index = frame_index; % frame_index by default
                    
                    % Use delayed index on uneven frames
                    if mod(j,2) == 1
                        index = delayed_frame_index;
                    end
  
                    output(index+symbol_index) = input(frame_index+symbol_index);
                        
                    % The middle 4 symbols up until 28 and last 4
                    % symbols need to be inverted
                    if (j >= 12 && j <= 15) || j>=28
                        output(index+symbol_index) = bitcmp(output(index+symbol_index),'uint8');
                    end    

                end
            end
            
            % As all symbols in all frames have now been delayed by 1
            % frame, the first frame is empty, so we remove it
            % Because the last frame's uneven symbols are empty and being
            % delayed into a new frame, the newly created new frame is also
            % empty and needs to be removed
            output = output(33:(n_frames)*32);
            n_frames = n_frames - 1;
            
        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_C1(obj,input,n_frames)
            % CIRC Decoder: C1 decoder
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            
            % The decoder will take away the 4 parity symbols leaving an
            % output of 28 symbols per frame
            output = zeros(n_frames*28,1,'uint8');
            
            % This will hold booleans for each symbol in all frames
            erasure_flags_out = zeros(n_frames*28,1,'logical');
            
            for i = 0:n_frames-1
                
                frame_index = i*32;
                new_frame_index = i*28;
                
                [dec_out,error] = step(obj.dec1,input(frame_index+1:frame_index+32), zeros(32,1,'logical'));
                
                % Adjust the flag for the symbol when erasure is detected
                % This is according to algorithm 1 in the assignment
                output(new_frame_index+1:new_frame_index+28) = dec_out;
                if ~(error == 0 || error == 1)
                    erasure_flags_out(new_frame_index+1:new_frame_index+28) = 1;
                end
            end
        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_delay_unequal(~,input,erasure_flags_in,n_frames)
            % CIRC Decoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder
            % -erasure_flags_out: the erasure flags at the output of this block (eraure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_unequal(obj.CIRC_enc_delay_unequal(x)) == x!
            
            % Similarly to the encoder, the symbols get delayed with
            % unequal delay. Now the first symbol gets the largest delay
            % whilst the last symbol gets no delay
            
            % The first symbol of the last frame will be delayed with 27
            % frames so we need to extend the output again.
            output = zeros((n_frames+27*4)*28,1,'uint8');
            % We need to make out erasures buffer correspond to the output
            erasure_flags_out = zeros((n_frames+27*4)*28,1,'logical');
            
            
            for i=0:n_frames-1
                
                frame_index = i * 28;
                
                % Traverse symbols in each frame in reverse order to gain
                % easy access to the corresponding delay
                for j=27:-1:0
                    
                    delay = j*4*28;
                    symbol_index = 28-j; %reverse order
                    
                    output(frame_index+delay+symbol_index) = input(frame_index+symbol_index);
                    erasure_flags_out(frame_index+delay+symbol_index) = erasure_flags_in(frame_index+symbol_index);
                end
            end
            
            % All of the symbols should now have been delayed with 27*4
            % frames which means the first 27*4 symbols are now empty and
            % can be removed
            output = output(27*4*28+1:(n_frames)*28);
            
            % Do the same operation on the erasure flags buffer
            erasure_flags_out = erasure_flags_out(27*4*28+1:(n_frames)*28);
            
            % Adjust the total number of frames
            n_frames = n_frames - 27*4;
            
        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_C2(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: C2 decoder
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            
            % prepare output and erasure flag buffers for use
            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
             
            
            
            for i=0:n_frames-1
                  
                frame_index = i*28;
                new_frame_index = i*24;

                % The parity symbols should now be in the middle of the
                % frame, to be able to decode we need to put them at the
                % end of the frame again.
                frame = zeros(28,1,'uint8');
                
                frame(1:12) = input(frame_index+1:frame_index+12);
                frame(13:24) = input(frame_index+17:frame_index+28);
                frame(25:28) = input(frame_index+13:frame_index+16);
                
                % Do the same for the erasure buffers frame
                erasure_frame = zeros(28,1,'logical');
                
                erasure_frame(1:12) = erasure_flags_in(frame_index+1:frame_index+12);
                erasure_frame(13:24) = erasure_flags_in(frame_index+17:frame_index+28);
                erasure_frame(25:28) = erasure_flags_in(frame_index+13:frame_index+16);
                
                % Perform the C2 decoding according to the algorithm
                % specified in the assignment
                % The amount of erasure flags at the input of c2
                erasure_count = sum(erasure_frame(:)==1);
                
                % Decoding can only be done when the amount of errors in a
                % frame passed from c1 is less or equal than 2
                if erasure_count <= 4
                    
                    % The decoder will strip the 4 parity bits and try to
                    % correct errors. The output frame will be 24 symbols in
                    % size
                    [dec_out, error] = step(obj.dec2,frame,erasure_frame);
                    
                    if error == 0 || error == 1
                        output(new_frame_index+1:new_frame_index+24) = dec_out;
                    else
                        if erasure_count > 2
                            erasure_flags_out(new_frame_index+1:new_frame_index+24) = erasure_frame(1:24);
                        
                        elseif erasure_count == 2 && error~=-1
                            output(new_frame_index+1:new_frame_index+24) = dec_out;
                        else 
                            erasure_flags_out(new_frame_index+1:new_frame_index+24) = 1;
                        end
                    end
                
                % in case there are more errors coming from c1 than the c2
                % can handle, just update all erasures for the frame
                else
                    % probably not needed but maybe better result
                    %output(new_frame_index+1:new_frame_index+24) = frame(1:24);
                    % Copy over the erasure flags from the input
                    erasure_flags_out(new_frame_index+1:new_frame_index+24) = 1;
                end
            end
        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_deinterleave_delay(~,input,erasure_flags_in,n_frames)
            % CIRC Decoder: De-interleaving sequence + delay of 2 frames
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder
            % -erasure_flags_out: the erasure flags at the output of this block (eraure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_deinterleave_delay(obj.CIRC_enc_delay_interleave(x)) == x!
            
            % Prepare buffers with 2 extra frames
            output = zeros((n_frames+2)*24,1,'uint8');
            erasure_flags_out = zeros((n_frames+2)*24,1,'logical');
            
            % add delay + do interleaving
            for i=0:n_frames-1

                frame_index = i*24;
                delayed_frame_index = (i+2)*24;
                
                frame = input(frame_index+1:frame_index+24);
                erasure_frame = erasure_flags_in(frame_index+1:frame_index+24);
                
                % Now the uneven words need to get delayed by 2
                
                output(frame_index+1:frame_index+2) = frame(1:2); %L6n+0
                erasure_flags_out(frame_index+1:frame_index+2) = erasure_frame(1:2);
                
                output(frame_index+3:frame_index+4) = frame(7:8); %R6n+0
                erasure_flags_out(frame_index+3:frame_index+4) = erasure_frame(7:8);
                
                output(delayed_frame_index+5:delayed_frame_index+6) = frame(13:14); %L6n+1
                erasure_flags_out(delayed_frame_index+5:delayed_frame_index+6) = erasure_frame(13:14);
                
                output(delayed_frame_index+7:delayed_frame_index+8) = frame(19:20); %R6n+1
                erasure_flags_out(delayed_frame_index+7:delayed_frame_index+8) = erasure_frame(19:20);
                
                output(frame_index+9:frame_index+10) = frame(3:4); %L6n+2
                erasure_flags_out(frame_index+9:frame_index+10) = erasure_frame(3:4);
                
                output(frame_index+11:frame_index+12) = frame(9:10); %R6n+2
                erasure_flags_out(frame_index+11:frame_index+12) = erasure_frame(9:10);
                
                
                
                output(delayed_frame_index+13:delayed_frame_index+14) = frame(15:16); %L6n+3
                erasure_flags_out(delayed_frame_index+13:delayed_frame_index+14) = erasure_frame(15:16);
                
                output(delayed_frame_index+15:delayed_frame_index+16) = frame(21:22); %R6n+3
                erasure_flags_out(delayed_frame_index+15:delayed_frame_index+16) = erasure_frame(21:22);
                
                output(frame_index+17:frame_index+18) = frame(5:6); %L6n+4
                erasure_flags_out(frame_index+17:frame_index+18) = erasure_frame(5:6);
                
                output(frame_index+19:frame_index+20) = frame(11:12); %R6n+4
                erasure_flags_out(frame_index+19:frame_index+20) = erasure_frame(11:12);
                
                output(delayed_frame_index+21:delayed_frame_index+22) = frame(17:18); %L6n+5
                erasure_flags_out(delayed_frame_index+21:delayed_frame_index+22) = erasure_frame(17:18);
                
                output(delayed_frame_index+23:delayed_frame_index+24) = frame(23:24); %R6n+5
                erasure_flags_out(delayed_frame_index+23:delayed_frame_index+24) = erasure_frame(23:24);
                
            end
            
            % Remove empty first and last 2 empty frames because they are
            % empty
            output = output(2*24+1 : n_frames*24);
            erasure_flags_out = erasure_flags_out(2*24+1 : n_frames*24);
            % Update the total number of frames
            n_frames = n_frames - 2;
            
        end
        
        function [output,n_frames] = C3_enc_8_parity(obj,input,n_frames)
            % Configuration 3: Generation of 8 parity symbols
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:n_frames
                output((i-1)*32+1:i*32) = step(obj.enc_8_parity,input((i-1)*24+1:i*24));
            end
        end
        
        function [output,erasure_flags_out,n_frames] = C3_dec_8_parity(obj,input,n_frames)
            % Configuration 3: Decoder
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block
            % -erasure_flags_out: the erasure flags at the output of this block
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            for i = 1:n_frames
                [output_dec,ERR] = step(obj.dec_8_parity,input((i-1)*32+1:i*32));
                if ERR == -1
                    output((i-1)*24+1:i*24) = output_dec;
                    erasure_flags_out((i-1)*24+1:i*24) = 1;
                else
                    output((i-1)*24+1:i*24) = output_dec;
                end
                
            end
        end
        
        function [output,interpolation_failed] = interpolator(obj,input,erasure_flags_in)
            % Interpolation: Linear interpolation
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block
            % -erasure_flags_in: the erasure flags at the input of this block
            % OUTPUT:
            % -output: linear interpolation of the input where there are no more than obj.max_interpolation consecutive erasures
            % -interpolation_failed: equal to one at the samples where interpolation failed
            
            % Make sure first and last sample are defined. If erasure: set to 0 (corresponds to the unsigned integer 2^15)
            if erasure_flags_in(1) ~= 0
                erasure_flags_in(1) = 0;
                input(1) = 2^15;
            end
            if erasure_flags_in(end) ~= 0
                erasure_flags_in(end) = 0;
                input(end) = 2^15;
            end
            
            erasure_burst = zeros(size(erasure_flags_in(:))); % Number of consecutive erasures
            ii = strfind([0,erasure_flags_in(:)'],[0 1]);
            erasure_burst(ii) = strfind([erasure_flags_in(:)',0],[1 0]) - ii + 1;
            
            output = input;
            interpolation_failed = erasure_flags_in;
            for i = find( (erasure_burst>0) & (erasure_burst<= obj.max_interpolation))'
                output(i:i+erasure_burst(i)-1) = uint16(round(double(output(i-1))+(1:erasure_burst(i))*(double(output(i+erasure_burst(i)))-double(output(i-1)))/(erasure_burst(i)+1)));
                interpolation_failed(i:i+erasure_burst(i)-1) = 0;
            end
            
            
        end
        
    end
    
    methods(Static)
        function test()
            % Test the Matlab code of this class
            
            audio_file = audioread('Hallelujah.wav');
            Fs = 44.1e3;
            
            %audio_file = audioread('Hallelujah_22050.wav');
            %Fs = 22.05e3;
            
            cd = AudioCD(Fs,1,8);
            cd = cd.writeCd(audio_file);
            
            T_scratch = 600000; % Scratch at a diameter of approx. 66 mm
            for i = 1:floor(numel(cd.cd_bits)/T_scratch)
                cd = cd.scratchCd(3000,30000+(i-1)*T_scratch);
            end
            cd = cd.bitErrorsCd(0.001);
            tic
            [out,interpolation_flags] = cd.readCd();
            toc
            
            fprintf('Number samples with erasure flags: %d\n',sum(sum(interpolation_flags~=0)));
            fprintf('Number samples with failed interpolations: %d\n',sum(sum(interpolation_flags==-1)));
            fprintf('Number undetected errors: %d\n',sum(sum(out(interpolation_flags==0) ~= cd.scaled_quantized_padded_original(interpolation_flags==0))));
            
            sound(out,Fs);
            
        end
        
        function test_scratchlengths()
           
            audio_file = audioread('Hallelujah.wav');
            Fs = 44.1e3;
            
            widths = [100, 3000, 10000];
            
            erasures = zeros(4,3);
            failed_interpol = zeros(4,3);
            undetected_errors = zeros(4,3);
            
            for configuration=0:3
                for i=1:3
                    
                    fprintf('config: %d, burst width: %d\n', configuration, widths(i));
                    
                    cd = AudioCD(Fs,configuration,8);
                    cd = cd.writeCd(audio_file);
                    % how many bits you need to jump over to get back to the same angle on the circular CD
                    T_scratch = 600000;
                    %amount of jumps
                    n_steps = floor( numel(cd.cd_bits)/T_scratch );
                    for j = 0:n_steps-1
                        cd = cd.scratchCd(widths(i),30000+(j)*T_scratch);
                    end
                    
                    [out,flags] = cd.readCd();
                    
                    %erasure, undetected error and failed interpolation probability for a
                    %certain config and a certain error probability
                    %flag is 0 for erasure and -1 for interpolation fail
                    erasures(configuration+1,i) = sum(sum(flags ~= 0)) / numel(flags);
                    failed_interpol(configuration+1,i) = sum(sum(flags == -1))/numel(flags);
                    undetected_errors(configuration+1,i) = sum(sum(out(flags==0) ~= cd.scaled_quantized_padded_original(flags==0)))/numel(flags);
                    
                end
            end
            
            erasures
            failed_interpol
            undetected_errors
            
        end
        
        
            
        
        function test_biterrorrate()
            
            audio_file = audioread('Hallelujah.wav');
            Fs = 44.1e3;
            %logspace(x, y, N) generates N samples logarithmic equally-spaced between
            %10**x and 10**y 
            errorlog = logspace(-1-log10(2),-3,10); %0.001 -> 0.05 
            
            %keep track of the flagged erasures
            erasures = zeros(3,10);
            failed_interpol = zeros(3, 10);
            for configuration=1:3
                for i=1:10
                    
                    fprintf('configuration %d\n',configuration);
                    fprintf('error probability %d\n',errorlog(i));
                    
                    cd = AudioCD(Fs,configuration,8);
                    cd = cd.writeCd(audio_file);
                    cd = cd.bitErrorsCd(errorlog(i));
                    
                    [out,flags] = cd.readCd();
                    
                    %erasure and failed interpolation probability for a
                    %certain config and a certain error probability
                    %flag is 0 for erasure and -1 for interpolation fail
                    erasures(configuration,i) = sum(sum(flags ~= 0)) / numel(flags);
                    failed_interpol(configuration,i) = sum(sum(flags == -1))/numel(flags);
                    
                end
            end
            errorlog
            erasures
            failed_interpol
            figure(1)
            semilogx(errorlog, erasures)
            grid on
            legend('config 1', 'config 2', 'config 3')
            figure(2)
            semilogx(errorlog, failed_interpol)
            grid on
            legend('config 1', 'config 2', 'config 3')
            
        end
        
    end
    
end
