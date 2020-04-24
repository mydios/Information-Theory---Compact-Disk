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
            
            % step 1: generating systematic code
            code = [msg, gf(zeros(size(msg,1), 2*obj.t), obj.m)];  
            
            for row=1:size(code,1)
                
                % step 2: generating systematic code
                [q, r] = deconv(code(row,:),obj.g);
                
                % step 3: generating systematic code
                code(row,:) = code(row,:)+r;
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
            
            % add zeros to the shortened code
            padding = [gf(zeros(size(code,1), obj.k-obj.l), obj.m), code];
            code = padding;
            
            decoded = gf(zeros(size(code,1), obj.l), obj.m);
            nERR = zeros(size(code,1),1);
            primel = gf(2,obj.m); % primitive element
            
            for row=1:size(code,1)
                
                % Get received codeword
                r = code(row, :);
                
                % Calculate syndrome polynomial (in order of descending powers)
                S = gf(zeros(1, 2*obj.t),obj.m);
                for syndrome=1:2*obj.t
                    S(syndrome) = S(syndrome) + polyval(r, primel^(2*obj.t-syndrome+obj.m0));
                end
                
                % allocate mem for code word we are looking for
                c = gf(zeros(1, obj.n),obj.m);
                
                % only proceed if the Syndrome polynomial has nonzero
                % values, otherwise r is a code word
                if any(S)
                    
                    % Do extended euclidean algorithm to obtain error locator
                    % polynomial

                    % step 1: initialize
                    Omega_2 = gf(zeros(1,2*obj.t+1),obj.m);
                    Omega_2(1) = gf(1,obj.m);
                    Omega_1 = S;
                    Delta_2 = gf(0,obj.m);
                    Delta_1 = gf(1,obj.m);

                    % step 2: iterate
                    stop = 0;
                    while stop == 0
                       [q, rest] = deconv(Omega_2, Omega_1);
                       Omega_2 = Omega_1;
                       
                       % get rid of zeros in the first powers
                       ind = find(rest~=0, 1, 'first');
                       Omega_1 = rest(ind:end);
                       Delta_temp = Delta_1;
                       
                       % make both arrays same size for addtition
                       mult = conv(q, Delta_1);
                       if size(mult,2) > size(Delta_2,2)
                           Delta_t = [gf(zeros(1,size(mult,2)-size(Delta_2,2)),obj.m), Delta_2];
                           Delta_2 = Delta_t;
                       elseif size(mult,2) < size(Delta_2,2)
                           mult_temp = [gf(1,zeros(size(Delta_2,2)-size(mult,2)),obj.m), mult];
                           mult = mult_temp;
                       end
                       Delta_1 = Delta_2 + mult;
                       Delta_2 = Delta_temp;
   
                       % step 3: terminate
                       if size(Omega_1, 2) < (obj.t+1)
                          stop = 1; 
                       end

                    end
                    
                    % step 4: make Delta have a constant term of 1
                    Delta = Delta_1;
                    constant = Delta(size(Delta,2));

                    exponent = 0;
                    if constant ~= 1
                       exponent = obj.n-log(constant);
                    end
                    Delta = conv(Delta, primel^(exponent));
                    
                    % Calculate Fourier transform E of the error vector e
                    v = size(Delta,2)-1;
                    E = [gf(zeros(1,obj.m0),obj.m),fliplr(S),gf(zeros(1,obj.n-obj.m0-2*obj.t),obj.m)];
                    if v > 0
                        for j=(obj.m0+2*obj.t+1):obj.n
                            mul = E(1,j-v:j-1).*Delta(1,1:v);
                            for el=1:v
                                E(1,j) = E(1,j) + mul(1,el);
                            end
                        end
                        for j=obj.m0:-1:1
                            temp = E(1,j+1:j+v);
                            mul = (Delta(1)^(-1))*(temp.*Delta(1,2:v+1));
                            for el=1:v
                                E(1,j) = E(1,j)+mul(1,el);
                            end
                        end
                    end

                    % Do inverse Fourier transform to obtain the error vector e
                    e = gf(zeros(1,obj.n),obj.m);
                    for i=0:(2*obj.t+obj.l-1)
                        resip = gf(zeros(1,obj.n),obj.m);
                        for j=0:(obj.n-1)
                            resip(1,j+1) = primel^(-j*i);
                        end
                        temp = E.*resip;
                        for j=1:(obj.n)
                            e(1,i+1) = e(1,i+1) + temp(1,j);
                        end
                    end
                    
                    % Determine the code word
                    c = r + fliplr(e);
                    
                    % Determine amount of errors
                    nErr_i = 0;
                    [div, rest] = deconv(c,obj.g);
                    if any(rest)
                        nErr_i = -1;
                    else
                        for i=1:(obj.n)
                           if e(1,i) ~= 0
                               nErr_i = nErr_i + 1;
                           end
                        end
                    end
                    nERR(row) = nErr_i;
                    
                else
                    % in case if all elements of S are zero, r is a code word
                    c = r;
                    
                    % no errors are found
                    nERR(row) = 0;
                end
               
                % Determine the information word
                % get rid of the zeros in the front and parity symbols in
                % the back
                b = c(1,1+obj.n-2*obj.t-obj.l:obj.n-2*obj.t);
                decoded(row,:) = b;
            end
        
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
                      
            primel = gf(2,m); % primitive element
            generator = [gf(1,m) primel^m0];
            for exponent = 1:(2*t-1)
                generator = conv(generator, [gf(1,m) primel^(m0+exponent)]);
            end
            
        end
        
        function h = determineCheckPolynomial(m, t, m0)
            % Generate the Reed-Solomon check polynomial with error correcting capability t over GF(2^m)
            % INPUT:
            % -m: positive integer
            % -t: error correction capability of the Reed-Solomon code, positive integer > 1
            % -m0: determines first root of the corresponging generator polynomial, positve integer >= 0
            % OUTPUT:
            % -h: rowvector containing the log of the GF(2^m) check polynomial coefficients in order of descending powers
            
            g = RSCode.makeGenerator(m, t, m0);
            x = gf(zeros(1,2^m),m);
            x(1) = gf(1,m);
            x(2^m) = gf(1,m);
            [h, r] = deconv(x,g);
            h = log(h);
            
        end
        
        function test()
            % Test the Matlab code of this class
            
            m0 = 1; % Also test with other values of m0!
            
            rs = RSCode(8,5,10,m0); % Construct the RSCode object
            
            msg = gf(randi([0,2^8-1],5,10),8); % Generate a random message of 5 information words
            
            code = rs.encode(msg); % Encode this message
            
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
        
        