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
                obj.dec2 = comm.RSDecoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);%'ErasuresInputPort',true);

                obj.gpoly1 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc1 = comm.RSEncoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec1 = comm.RSDecoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);%'ErasuresInputPort',true);
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
    
        function [output,n_frames] = CIRC_enc_delay_interleave(obj,input,n_frames)
            % CIRC Encoder: Delay of 2 frames + interleaving sequence
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % 2 delay frames are added to the output
            delay = 2;
            n_frames = n_frames + delay;
            output = zeros(n_frames*24,1,'uint8');
            
            % add delay + do interleaving
            for i=1:(n_frames-delay)
                
                frame = input((i-1)*24+1:i*24);
                
                output((i+1)*24+1:(i+1)*24+2) = frame(1:2);
                output((i+1)*24+7:(i+1)*24+8) = frame(3:4);
                output((i+1)*24+3:(i+1)*24+4) = frame(9:10);
                output((i+1)*24+9:(i+1)*24+10) = frame(11:12);
                output((i+1)*24+5:(i+1)*24+6) = frame(17:18);
                output((i+1)*24+11:(i+1)*24+12) = frame(19:20);
                
                output((i-1)*24+13:(i-1)*24+14) = frame(5:6);
                output((i-1)*24+19:(i-1)*24+20) = frame(7:8);
                output((i-1)*24+15:(i-1)*24+16) = frame(13:14);
                output((i-1)*24+21:(i-1)*24+22) = frame(15:16);
                output((i-1)*24+17:(i-1)*24+18) = frame(21:22);
                output((i-1)*24+23:(i-1)*24+24) = frame(23:24);
                
            end
            
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
            
            output = zeros(n_frames*28,1,'uint8');
            for i = 1:n_frames
                encoded = step(obj.enc2,input((i-1)*24+1:i*24));
                output((i-1)*28+1:(i-1)*28+12) = encoded(1:12);
                output((i-1)*28+13:(i-1)*28+16) = encoded(25:28);
                output((i-1)*28+17:i*28) = encoded(13:24);
            end
            
        end
        
        function [output,n_frames] = CIRC_enc_delay_unequal(obj,input,n_frames)
            % CIRC Encoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            D = 4;
            output = zeros((n_frames+27*D)*28,1,'uint8');
            
            for i=1:n_frames
               for j=0:27
                  output((i-1+j*D)*28+j+1) = input((i-1)*28+j+1); 
               end
            end
            n_frames = n_frames + 27*D;
            
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
            
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:n_frames
                output((i-1)*32+1:i*32) = step(obj.enc1,input((i-1)*28+1:i*28));
            end
            
        end
        
        function [output,n_frames] = CIRC_enc_delay_inv(obj,input,n_frames)
            % CIRC Encoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            output = zeros((n_frames+1)*32,1,'uint8');
            for i=1:n_frames
                for j=1:32
                   
                   if mod(j,2) == 1
                       if (j >= 13 && j <= 16) || j>=29
                           output(i*32+j) = bitcmp(input((i-1)*32+j),'uint8');
                       else
                           output(i*32+j) = input((i-1)*32+j);
                       end
                   else
                       if (j >= 13 && j <= 16) || j>=29
                           output((i-1)*32+j) = bitcmp(input((i-1)*32+j),'uint8');
                       else
                           output((i-1)*32+j) = input((i-1)*32+j);
                       end
                   end 
                end
            end
            n_frames = n_frames + 1;
            
        end
        
        function [output,n_frames] = CIRC_dec_delay_inv(obj,input,n_frames)
            % CIRC Decoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_inv(obj.CIRC_enc_delay_inv(x)) == x!

            output = zeros((n_frames+1)*32,1,'uint8');
            for i=1:n_frames
                for j=1:32
                   
                   if mod(j,2) == 0
                       if (j >= 13 && j <= 16) || j>=29
                           output(i*32+j) = bitcmp(input((i-1)*32+j),'uint8');
                       else
                           output(i*32+j) = input((i-1)*32+j);
                       end
                   else
                       if (j >= 13 && j <= 16) || j>=29
                           output((i-1)*32+j) = bitcmp(input((i-1)*32+j),'uint8');
                       else
                           output((i-1)*32+j) = input((i-1)*32+j);
                       end
                   end 
                end
            end
            
            % get rid of first and last frame (they are empty)
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

            output = zeros(n_frames*28,1,'uint8');
            erasure_flags_out = zeros(n_frames*28,1,'logical');
            
            for i = 1:n_frames
                
                [output_dec,ERR] = step(obj.dec1,input((i-1)*32+1:i*32));%,zeros(32,1,'logical'));
                if ERR == 0 || ERR == 1
                    output((i-1)*28+1:i*28) = output_dec;
                else
                    output((i-1)*28+1:i*28) = output_dec;
                    erasure_flags_out((i-1)*28+1:i*28) = 1;
                end
            end
            
        end

        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_delay_unequal(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (erasure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_unequal(obj.CIRC_enc_delay_unequal(x)) == x!
            
            D = 4;
            output = zeros((n_frames+27*D)*28,1,'uint8');
            erasure_flags_out = zeros((n_frames+27*D)*28,1,'logical');
            
            for i=1:n_frames
               for j=27:-1:0
                  output((i-1+j*D)*28+28-j) = input((i-1)*28+28-j);
                  erasure_flags_out((i-1+j*D)*28+28-j) = erasure_flags_in((i-1)*28+28-j);
               end
            end
            
            % get rid of first and last 27*D frames (they are empty)
            output = output(27*D*28+1:(n_frames)*28);
            erasure_flags_out = erasure_flags_out(27*D*28+1:(n_frames)*28);
            n_frames = n_frames - 27*D;
            
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
            
            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            
            for i=1:n_frames
               
                % put parity symbols in the back of the frame
                encoded = zeros(28,1,'uint8');
                encoded(1:12) = input((i-1)*28+1:(i-1)*28+12);
                encoded(13:24) = input((i-1)*28+17:i*28);
                encoded(25:28) = input((i-1)*28+13:(i-1)*28+16);
                
                erasure_encoded = zeros(28,1,'logical');
                erasure_encoded(1:12) = erasure_flags_in((i-1)*28+1:(i-1)*28+12);
                erasure_encoded(13:24) = erasure_flags_in((i-1)*28+17:i*28);
                erasure_encoded(25:28) = erasure_flags_in((i-1)*28+13:(i-1)*28+16);
                
                % decode
                [decoded, ERR] = step(obj.dec2,encoded);%,erasure_encoded);
                
                % amount of erasure flags at the input
                f = sum(erasure_encoded(:)==1);
                
                % algorithm for C2 decoding
                if ERR == 0 || ERR == 1
                    output((i-1)*24+1:i*24) = decoded;
                else
                   if  f > 2
                       output((i-1)*24+1:i*24) = encoded(1:24);
                       erasure_flags_out((i-1)*24+1:i*24) = erasure_encoded(1:24);
                   elseif f == 2 && ERR ~= -1
                       output((i-1)*24+1:i*24) = decoded;
                   else
                       output((i-1)*24+1:i*24) = decoded;
                       erasure_flags_out((i-1)*24+1:i*24) = 1;
                   end
                end
                
            end
           
        end
        
        
        function [output,erasure_flags_out,n_frames] = CIRC_dec_deinterleave_delay(obj,input,erasure_flags_in,n_frames)
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
            
            delay = 2;
            output = zeros((n_frames+delay)*24,1,'uint8');
            erasure_flags_out = zeros((n_frames+delay)*24,1,'logical');
            
            % add delay + do interleaving
            for i=1:n_frames
               
               frame = input((i-1)*24+1:i*24);
               erasures = erasure_flags_in((i-1)*24+1:i*24);
               
               output((i-1)*24+1:(i-1)*24+2) = frame(1:2);
               erasure_flags_out((i-1)*24+1:(i-1)*24+2) = erasures(1:2);
               output((i-1)*24+3:(i-1)*24+4) = frame(7:8);
               erasure_flags_out((i-1)*24+3:(i-1)*24+4) = erasures(7:8);
               output((i-1)*24+9:(i-1)*24+10) = frame(3:4);
               erasure_flags_out((i-1)*24+9:(i-1)*24+10) = erasures(3:4);
               output((i-1)*24+11:(i-1)*24+12) = frame(9:10);
               erasure_flags_out((i-1)*24+11:(i-1)*24+12) = erasures(9:10);
               output((i-1)*24+17:(i-1)*24+18) = frame(5:6);
               erasure_flags_out((i-1)*24+17:(i-1)*24+18) = erasures(5:6);
               output((i-1)*24+19:(i-1)*24+20) = frame(11:12);
               erasure_flags_out((i-1)*24+19:(i-1)*24+20) = erasures(11:12);
               
               output((i+1)*24+5:(i+1)*24+6) = frame(13:14);
               erasure_flags_out((i+1)*24+5:(i+1)*24+6) = erasures(13:14);
               output((i+1)*24+7:(i+1)*24+8) = frame(19:20);
               erasure_flags_out((i+1)*24+7:(i+1)*24+8) = erasures(19:20);
               output((i+1)*24+13:(i+1)*24+14) = frame(15:16);
               erasure_flags_out((i+1)*24+13:(i+1)*24+14) = erasures(15:16);
               output((i+1)*24+15:(i+1)*24+16) = frame(21:22);
               erasure_flags_out((i+1)*24+15:(i+1)*24+16) = erasures(21:22);
               output((i+1)*24+21:(i+1)*24+22) = frame(17:18);
               erasure_flags_out((i+1)*24+21:(i+1)*24+22) = erasures(17:18);
               output((i+1)*24+23:(i+1)*24+24) = frame(23:24);
               erasure_flags_out((i+1)*24+23:(i+1)*24+24) = erasures(23:24);
                
            end
            
            % git rid of first and last 2 frames (they are empty)
            output = output(delay*24+1 : n_frames*24);
            erasure_flags_out = erasure_flags_out(delay*24+1 : n_frames*24);
            n_frames = n_frames - delay;
            
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
            %Test the Matlab code of this class
            
%             audio_file = audioread('Hallelujah.wav');
%             Fs = 44.1e3;
            
            audio_file = audioread('Hallelujah.wav');
            %Fs = 22.05e3;
            Fs = 44.1e3;


            cd = AudioCD(Fs,1,8);
            tic
            cd = cd.writeCd(audio_file);
            toc
            
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
        
        function Q5_1()
           
            audio_file = audioread('Hallelujah.wav');
            Fs = 44.1e3;
            
            burst = [100, 3000, 10000];
            
            n_erasures = zeros(4,3);
            n_failed = zeros(4,3);
            n_undetected = zeros(4,3);
            
            for configuration=0:3
                for i=1:3
                    
                    fprintf('configuration %d\n',configuration);
                    fprintf('burstlength %d\n',burst(i));
                    
                    cd = AudioCD(Fs,configuration,8);
                    tic
                    cd = cd.writeCd(audio_file);
                    toc
                    
                    T_scratch = 600000; % Scratch at a diameter of approx. 66 mm
                    for j = 1:floor(numel(cd.cd_bits)/T_scratch)
                        cd = cd.scratchCd(burst(i),30000+(j-1)*T_scratch);
                    end
                    
                    tic
                    [out,interpolation_flags] = cd.readCd();
                    toc
                    
                    n_erasures(configuration+1,i) = sum(sum(interpolation_flags~=0));
                    n_failed(configuration+1,i) = sum(sum(interpolation_flags==-1));
                    n_undetected(configuration+1,i) = sum(sum(out(interpolation_flags==0) ~= cd.scaled_quantized_padded_original(interpolation_flags==0)));
                    
                end
            end
            
            n_erasures
            n_failed
            n_undetected
            
        end
        
        function Q5_2()
            
            audio_file = audioread('Hallelujah.wav');
            Fs = 44.1e3;
            
            errorp = logspace(-1-log10(2),-3,10);
            
            n_erasures = zeros(3,10);
            n_failed = zeros(3,10);
            
            for configuration=1:3
                for i=1:10
                    
                    fprintf('configuration %d\n',configuration);
                    fprintf('error probability %d\n',errorp(i));
                    
                    cd = AudioCD(Fs,configuration,8);
                    tic
                    cd = cd.writeCd(audio_file);
                    toc
                    
                    cd = cd.bitErrorsCd(errorp(i));
                    
                    tic
                    [out,interpolation_flags] = cd.readCd();
                    toc
                    
                    n_erasures(configuration,i) = sum(sum(interpolation_flags~=0))/numel(interpolation_flags);
                    n_failed(configuration,i) = sum(sum(interpolation_flags==-1))/numel(interpolation_flags);
                    
                end
            end
            
            n_erasures
            n_failed
            
        end
        
    end
    
end