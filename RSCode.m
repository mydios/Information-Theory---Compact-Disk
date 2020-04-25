classdef RSCode
    % RSCode class: allows encoding and decoding using a Reed-Solomon code
    %
    % Author: Johannes Van Wonterghem, Jan 2017
    %
    % See the static test() method for an example usage of this class
    
    properties
        m; % GF(2^m) field
        
        n; % Code length
        k; % Information length
        t; % Error correction capability
        l; % Shortened information length (-> shortened code length = l+n-k)
        
        m0; % m0 of the Reed-Solomon code, determines first root of generator
        
        g; % rowvector containing the GF(2^m) generator polynomial coefficients in order of descending powers
    end
    
    methods
    
        function obj=RSCode(m,t,l,m0)
			% Constructor of the RSCode class
			% INPUT:
            % -m: Defines the GF(2^m) field
            % -t: Error correction capability
			% -l: Shortened information length (-> shortened code length = l+n-k)
			% -m0: m0 of the Reed-Solomon code, determines first root of generator
            % OUTPUT:
            % -obj: the RSCode object
            obj.m = m;
          
            obj.t = t;
            obj.l = l;
            
            obj.n = 2^m-1;
            obj.k = obj.n-2*t;
           
            obj.m0 = m0;
            
            obj.g = RSCode.makeGenerator(m,t,m0);

        end
        
        function code = encode(obj,msg)
            % Systematically encode information words using the Reed-Solomon code
            % INPUT:
            % -obj: the RSCode object, defines the properties of the code
            % -msg: every row contains a length obj.l information word consisting of GF(2^m) elements
            % OUTPUT:
            % -code: every row contains a GF(2^m) codeword corresponding to systematic Reed-Solomon coding of the corresponding information word            
            
            assert(size(msg,2) == obj.l);
            code = gf(zeros(size(msg,1),obj.l+obj.n-obj.k),obj.m);
            
            %encode word for word using second method on pg 145
            for i = 1:size(msg,1)
                msg_ex = gf([msg(i,:) zeros(1,obj.n-obj.k)],obj.m); %*x^(n-k)
                [~,rem] = deconv(msg_ex, obj.g); %divide by g(x)
                code(i,:) = msg_ex + rem; %add remainder to dividend to obtain codeword
            end
            
                        
        end
        
        function [decoded,nERR] = decode(obj,code)
            % Decode Reed-Solomon codes
            % INPUT:
            % -obj: the RSCode object, defines the properties of the code
            % -code:  every row contains a GF(2^m) codeword of length obj.l+2*obj.t
            % OUTPUT:
            % -decoded: every row contains a GF(2^m) information word corresponding to decoding of the corresponding Reed-Solomon codeword
            % -nERR: column vector containing the number of corrected symbols for every codeword, -1 if error correction failed
            
            assert(size(code,2) == obj.l+2*obj.t);
            nERR = zeros(size(code, 1),1);
            ext_code = gf(zeros(size(code,1), obj.n),obj.m);
            corr_code = gf(zeros(size(code,1), obj.n),obj.m);
            four_code = gf(zeros(size(code, 1), obj.n), obj.m);
            ext_code(:,(obj.n-obj.l-2*obj.t)+1:obj.n) = code;
            temp = ext_code.x;
            alpha = gf(2,obj.m);
            
            
            aij = dftmtx(alpha);
            aijinv = dftmtx(1/alpha);
            for i = 1:size(code,1)
                
                %find syndrome polynome
                S = polyval(gf(temp(i,:),obj.m),alpha.^(2*obj.t-1+obj.m0:-1:obj.m0));                
                if(S == gf(zeros(1,2*obj.t),obj.m))
                    nERR(i) = 0;
                    corr_code(i,:) = ext_code(i,:);
                    continue;
                end
                %find Lambda
                prevLambda = gf(0,obj.m);
                Lambda = gf(1, obj.m);
                prevOmega = gf(zeros(1, 2*obj.t+1), obj.m);
                prevOmega(1) = gf(1, obj.m);
                Omega=S;
                while(size(Omega,2)-1 >= obj.t)
                    [q, nextOmega] = deconv(prevOmega, Omega);
                    c = conv(q, Lambda); 
                    nextLambda = [gf(zeros(1,max(length(c)-length(prevLambda),0)), obj.m) prevLambda] - [gf(zeros(1,max(length(prevLambda)-length(c),0)), obj.m) c];
                    prevLambda = Lambda;
                    Lambda = nextLambda(find(nextLambda ~= 0,1,'first'):end);
                    prevOmega = Omega;
                    Omega = nextOmega(find(nextOmega ~= 0,1,'first'):end);
                end
                %normalize
                Omega = Omega/Lambda(end);
                Lambda = Lambda/Lambda(end);
                Lambda = Lambda(find(Lambda ~= 0,1,'first'):end);
                if size(Lambda,2)-1 <= obj.t 
                    nERR(i) = size(Lambda,2)-1;
                else
                    nERR(i) = -1;
                    corr_code(i,:) = ext_code(i,:);
                    continue;
                end
                
                %find error vector in fourier domain
                E = [gf(zeros(1,obj.n-obj.m0-2*obj.t),obj.m) S gf(zeros(1,obj.m0),obj.m)];
                nu = size(Lambda,2)-1;
                for j = obj.n-obj.m0+1:obj.n
                    E(j) = E(j-nu:j-1)*Lambda(end:-1:2).'/Lambda(1); 
                end
                for j = obj.n-obj.m0-length(S):-1:1
                    E(j) = E(j+1:j+nu)*Lambda(end-1:-1:1).';                    
                end
                      
                %inverse fourier transform
                e = fliplr(E)*aijinv;
                e = fliplr(e);
                
                %correct codewords with e_hat
                if (e(1:end-obj.l-obj.n+obj.k) == gf(zeros(1,obj.k-obj.l),obj.m))
                    corr_code(i,:) = ext_code(i,:) + e;
                else
                    corr_code(i,:) = ext_code(i,:);
                    nERR(i) = -1;
                end
                
            end
            decoded = corr_code(:,end-obj.l-2*obj.t+1:end-2*obj.t);
        end
    
    end
    
    methods(Static)
        
        function generator = makeGenerator(m, t, m0)
            % Generate the Reed-Solomon generator polynomial with error correcting capability t over GF(2^m)
            % INPUT:
            % -m: positive integer
            % -t: error correction capability of the Reed-Solomon code, positive integer > 1
            % -m0: determines first root of generator polynomial, positve integer >= 0
            % OUTPUT:
            % -generator: rowvector containing the GF(2^m) generator polynomial coefficients in order of descending powers
            % generator = rsgenpoly(7,3);
            alpha = gf(2,m);
            gengf = 1;
            for k = m0:(m0+2*t-1)
                gengf = conv(gengf, [1 -(alpha.^k)]);
            end
            generator = gengf.x;
            
            
            
        end
        
        function test()
            % Test the Matlab code of this class
            
            m0 = 1; % Also test with other values of m0!
            
            rs = RSCode(8,5,10,m0); % Construct the RSCode object
            
            msg = gf(randi([0,2^8-1],5,10),8); % Generate a random message of 5 information words
            
            code = rs.encode(msg); % Encode this message
            tempcode = code;
            % Introduce errors
            code(2,[3 18]) = code(2,[5 18])+1;
            code(3,8) = 0;
            code(4,[4 2 19 20 6]) = gf(randi([0,2^8-1],1,5),8);
            code(5,[4 2 19 20 6 13]) = gf(randi([0,2^8-1],1,6),8);
            
            
            [decoded,nERR] = rs.decode(code); % Decode
            
            nERR
            
            assert(all(all(decoded(1:4,:) == msg(1:4,:))))
            
        end
        
    end
    
    
end
        
        