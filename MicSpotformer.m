%Author:    Dimme de Groot
%Date:      Sept. 2024
%Descr:     This object implements the microphone spotformer
%
%Methods:   
%       Create the spotformer object:           
%           MicSpot = MicSpotformer(c, fs, window_length, pad_length, N_int, IntWinRad, TarWinRad, nSigma2, numSigma2, rebRatio, flag_full_axis, analysis_window, synthesis_window;   
%               with:   c                   speed of sound [m/s], 
%                       fs                  sampling frequency [Hz]
%                       window_length       the length of the input frame in [s], 
%                       pad_length          the length of the padding in [s],
%                       N_int               the number of points per integral (triple integral: N_int^3 points in total) 
%                       IntWinRad           the 3 sigma radius of the interfer window
%                       TarWinRad           the 3 sigma radius of the target window
%                       nSigma2             a regularisation term which can be used in case of micropone self noise
%                       numSigma2           a regularisation term which is used to ensure the matrices are positive definite in the presence of numerical inaccuracies
%                       rebRatio            a term specifying the direct to reverberant sound ratio
%                       flag_full_axis      specifies whether or not to use the full frequency axis or only half (true/false)
%                       analysis_window     a string specifying the analysis window (currently only: "sqrthann") 
%                       synthesis_window    a string specifying the synthesis window (currently only: "sqrthann") 
%               Unless otherwise mentioned, all inputs are scalars!!!          
%           
%           See example files for the other methods    
%   
%Dependencies:  The function "clenquad" is needed. This function computes the quadrature weights. I used (a renamed version of) the following function:  
%                   https://www.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature
%               It has a permissive license and is provided on the github
%               
%               Note that the RIR-generator by E. Habets is provided on the github as well. This is not a dependency for the spotformer object itself, but is needed to run the examples.

classdef MicSpotformer < handle
    properties
        c                   %[m/s], speed of sound
        fs                  %[Hz],  sampling frequency

        window_length_act   %[s], actual window length
        N_t                 %[-], window length in samples      
        pad_length_act      %[s], actual padding length
        N_pad               %[-], padding length in samples
        N_hop               %[-], hop length in samples        
        w_analysis          %[-], analysis window
        w_synthesis         %[-], synthesis window  

        N_int;

        numSigma2           %[], regularisation term for numerical inaccuracy. In my experience, this can be zero. But very small will enforce positive definiteness in the presence of numerical inaccuracies
        nSigma2             %[], term for microphone self-noise 
        rebRatio            %[-], Parameter of regularisation applied to estimated covariance matrics
       
        Riso                %[], the covariance matrix corresponding to the late reverberation
        Rnum                %[], a regularisation term for mnumerical inaccuracies
        Rn                  %[], a regularisation term for self-noise
            
        R_Int               %[], the covariance matrix for the interferers
        R_Tar               %[], the covariance matrix for the target location 

        flag_full_axis      %[-], True for full frequency axis [-Fs/2, Fs/2). False for [0, Fs/2]
        k_ax                %[rad/m], wavenumber axis
        N_k                 %[-], length of k_axis

        IntWinSigma2        %[m2], sigma^2, interferer region
        TarWinSigma2        %[m2], sigma^2, target region

        weights_mic         %[-], the weights of the microphone spotformer
    end

    methods
        function obj = MicSpotformer(c, fs, window_length, pad_length, N_int, IntWinRad, TarWinRad, nSigma2, numSigma2, rebRatio, flag_full_axis, analysis_window, synthesis_window)
            if nargin == 0
                c = 342;
                fs = 16000;

                window_length = 0.016;
                pad_length = 0.016;
                    
                IntWinRad = 0.5;
                TarWinRad = 0.5;
                
                numSigma2 = 1e-10;
                nSigma2 = 0;
                
                rebRatio = 0;

                flag_full_axis = false;      
                analysis_window = "sqrthann";
                synthesis_window = "sqrthann";
            end
            obj.c = c;
            obj.fs = fs;
            obj.N_int = N_int;

            obj.N_t = 2^nextpow2(floor(fs.*window_length));    
            if pad_length ~= 0
                obj.N_pad = 2^nextpow2(floor(fs.*pad_length));   
            else
                obj.N_pad = 0;
            end
            obj.window_length_act = obj.N_t/fs;
            obj.pad_length_act = obj.N_pad/fs;
            
            obj.IntWinSigma2 = (IntWinRad/3)^2;         
            obj.TarWinSigma2 = (TarWinRad/3)^2;      

            obj.rebRatio = rebRatio;

            obj.nSigma2 = nSigma2;
            obj.numSigma2 = numSigma2;
            
            obj.flag_full_axis = flag_full_axis;    
            
            if analysis_window == "sqrthann"
               obj.N_hop = obj.N_t/2;
               obj.w_analysis = methodSqrthann(obj); 
            else
                disp("currently only sqrthann is supported as window :S")
            end

            if synthesis_window == "sqrthann"
               obj.w_synthesis = methodSqrthann(obj); 
            else
                disp("currently only sqrthann is supported as window :S")
            end

            methodK_ax(obj);
            obj.N_k = length(obj.k_ax);
        end

        %Some methods used internally

        function window = methodSqrthann(obj)
            % This function computes thes sqaure root hanning window
            % See, i.e., Smith, J.O. Spectral Audio Signal Processing, http://ccrma.stanford.edu/~jos/sasp/, online book, 2011 edition.
            n = (0:1:obj.N_t-1).'; %note: obj.N_t is the length of the window
            m = n-obj.N_t/2;
            window = sqrt(ones(obj.N_t, 1).*(0.5+0.5*cos(2*pi/obj.N_t * m)));
        end

        function methodK_ax(obj)
            % this function computes the wavenumber axis
            obj.N_k = obj.N_t + obj.N_pad;                              %[-], total number of frequency bins (even)  
            if obj.flag_full_axis
                k_ax = (-obj.N_k/2:1:obj.N_k/2-1)/obj.N_k*2*pi*obj.fs/obj.c; 
                k_ax = fftshift(k_ax); 
            else
                k_ax = (0:1:obj.N_k/2)/obj.N_k*2*pi*obj.fs/obj.c;          %[rad/m], wave number -> only half of the axis [0, Fs/2] is relevant
            end
            obj.k_ax = k_ax;
        end

        function fnc_comp_Rn(obj, dim)
            % This function computes Rnum and Rn
            obj.Rnum = eye(dim)*obj.numSigma2;
            obj.Rn = eye(dim)*obj.nSigma2;
        end

        function fnc_comp_Risotropic(obj, x)
            % This function computes Riso, excluding(!) the rebratio term.
            N = size(x,1);          %number of positions considered
            
            R_iso = zeros(N, N, length(obj.k_ax));
            dist = distcalc(x,x);   %norm(x1-x2), norm(x1-x2), etc.
        
            for i = 1:length(obj.k_ax)
                R_iso(:,:,i) = sinc(obj.k_ax(i)*dist);
            end
            obj.Riso = R_iso;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These are methods relating to the spatial integration %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function D=distcalc(X1,X2)
            %Descr: D = distcalc(X1, X2) computes the distances between X1 and X2 for each item in the list
            %       X1 is an M1 x dim matrix with entries (x,y,z) or similar (dim=3 for {x,y,z}, dim=2 for {x,y}, etc.)
            %       X2 is an M2 x dim "                                                                   " 
            %       D is an M1 x M2 matrix with the distances between the locations in X1 and X2. The i'th row of the j'th columnt corresponds to L2-norm( {X1}_i - {X2}_j ) 
            M1 = size(X1,1);
            M2 = size(X2,1);
            D = zeros(M1,M2);
            for i=1:M1
                D(i,:) = vecnorm(X1(i,:) - X2,2,2);
            end
        end

        function out = fnc_green_wavefunction(x_s, x_r, k)    
            %Descr: This function evaluates the greens function solution to the (acoustic) wave equation in the frequency domain.
            %           g(x_s, x_r, k) = exp(-1j * k * norm(x_s-x_r)) / (4*pi*norm(x_s-x_r)), with k = f/c the wavenumber, f the frequency in hertz, c the wave velocity
            %       Inputs: 
            %           x_s a 3x1 real vector giving the source location
            %           x_r a 3xN real vector giving the receiver location
            %           k   a real scalar giving the wavenumber
            %       Outputs:
            %           out a 1 x N complex vector. The i'th element contains the greens function evaluated for the i'th coordiante of x_r (i.e. x_r(:,i));
            %           
            %       Note: the greens function for the wave equation (actually: helmholtz equation) is symmetric in the coordinate argument. 
            %             I.e.  g(x_s, x_r, k) = g(x_r, x_s, k)

            dist = vecnorm(x_s-x_r);
            out = 1./(4*pi*dist).*exp(-1j*k*dist);
        end

        function out = fnc_integrand_Gaussian(x_s1, x_s2, x_r, k, mu_r, sigma_x, sigma_y, sigma_z, flag_spatial_weight)
            %Descr: This function is used for computing the spatial covariance matrices. 
            %           out = fnc_integrand_Gaussian(x_s1, x_s2, x_r, k, mu_r, sigma_x, sigma_y, sigma_z, flag_spatial_weight)
            %       If flag_spatial_weight==false     
            %           out = 1/(sqrt{(2*pi)^3}*sigma_x*sigma_y*sigma_z) * g(x_s1, x_r, k)*conj(g(x_s2, x_r, k))
            %       If flag_spatial_weight==true, a Gaussian is added:
            %           out = 1/(sqrt{(2*pi)^3}*sigma_x*sigma_y*sigma_z) * g(x_s1, x_r, k)*conj(g(x_s2, x_r, k))*exp(-0.5 {[(x_r(1)-mu(1))/sigma_x]^2+[(x_r(2)-mu(2))/sigma_y]^2+[(x_r(3)-mu(3))/sigma_z]^2})
            %
            %       Here, g is the greens function solution to the wave equation
            %                      
            %           Inputs: 
            %               - x_s1, x_s2 the source locations. 3x1 real vectors
            %               - x_r the receiver location. (integration variable) 3xN real vector --> the output is computed for each N
            %               - mu_r the mean of the Guassian. 3x1 real vector
            %               - sigma_x, sigma_y, sigma_z the standard deviations. Real scalar
            %               - k the wave number. Real scalar
            %               - flag_spatial_weight selects if the Gaussian weighting should be added.
            %           Output:
            %               - Out: 1xN complex vector. Each element out(1,i) is the result evaluated for receiver location x_r(:,i);

            %The inner term of the integrand
            out = fnc_green_wavefunction(x_s1, x_r, k).*conj(fnc_green_wavefunction(x_s2, x_r, k));
            
            %Normalisation term of (Gaussian) weighting
            normConst = sqrt((2*pi)^3)*sigma_x*sigma_y*sigma_z;
            
            %Gaussian (spatial) weighting
            %Note: if a Gauss-hermite quadrature is used, this is not needed (it is implicit in the weights). For Clenshaw-Curtis it is not included in the weights.  
            if flag_spatial_weight     
                diff = x_r - mu_r;
                out = out.*exp(-0.5*( (diff(1,:)/sigma_x).^2 + (diff(2,:)/sigma_y).^2  + (diff(3,:)/sigma_z).^2 ));
            end
        
            %Output should be normalised by normalisation term
            out = out/normConst;
        end

        function R_region = fnc_comp_Rregion(obj, x, x_bar, sigma2)
            %Descr:     Computes the covariance matrices over the regions. (Excluding the isotropic and numerical covariance matrices)     
            %
            %Inputs:    - N:        [1 x 1] or [3 x 1] or [1 x 3] the number of points per dimension in the numerical integration. Can easily be changed to make number of points dimension depedent
            %           - k_ax:     [Nk x 1], the frequency axis in radians per second 
            %           - x:        [Nr x 3], the locations for which the correlations are computed. 
            %           - x_bar:    [1 x 3],  the mean location around which we integrate (i.e. the integration variable, but we set the variable part in this function)
            %           - sigma_2:  [1 x 1] or [3 x 1] or [1 x 3], the standard deviation of the sphere we integrate. Can easily be changed to a vevctor for each dimension.  
            %Output:    - R_region: [Nr x Nr x Nk], the computed covariance matrix per frequency bin.  
            %
            %Dependencies: clenquad.m: a matlab function for computing the clenshaw-curtis quadrature weights.

            % Assign the number of integration points and standard deviations for the microphone spotformer
            
            N = obj.N_int;
            if isscalar(N)
                Nx = N; Ny = N; Nz = N;
            else
                Nx = N(1); Ny = N(2); Nz = N(3); 
            end
        
            % Assign the standard deviations per argument
            sigma = sqrt(sigma2);
            if isscalar(sigma)
                sigma_x = sigma; sigma_y = sigma; sigma_z = sigma; 
            end
        
            Nr = size(x,1);     %The number of microphones
            Nk = length(obj.k_ax);  %The number of frequency bins
        
            % Compute quadrature points: integrate from the mean location to +- 3 standard deviations
            [x_quad, w_x] = clenquad(Nx, -3*sigma_x+x_bar(1), 3*sigma_x+x_bar(1));
            [y_quad, w_y] = clenquad(Ny, -3*sigma_y+x_bar(2), 3*sigma_y+x_bar(2)); 
            [z_quad, w_z] = clenquad(Nz, -3*sigma_z+x_bar(3), 3*sigma_z+x_bar(3));
            [weight_list, coor_list] = weightlist(w_x, x_quad, w_y, y_quad, w_z, z_quad);
        
            % Lots of for loops! x-y-z dimension + frequency bins 
            % Note: R(i,j,k) = conj(R(j,i,k)): this allows for speed up!
            R_region = zeros(Nr, Nr, Nk);
            for i=1:Nr
                for j=1:Nr
                    for k=1:Nk
                        kk = obj.k_ax(k);  %Get frequency in rad/s        
        
                        out = fnc_integrand_Gaussian(x(i,:)', x(j,:)', coor_list, kk, x_bar', sigma_x, sigma_y, sigma_z, true);
                        R_region(i,j,k) = sum(weight_list.*out);   
                    end
                end
            end
        end

        function [weight_list, coor_list] = weightlist(w_x, coor_x, w_y, coor_y, w_z, coor_z)
            % Function to gert all coordinates and the corresponding weights needed for the numerical integration     
            %   Inputs:     - w_{x,y,z}:    weightvectors of length N_{x,y,z}. Orientation does not matter
            %               - coor_{x,y,z}: the corresponding coordinates of length  N_{x,y,z}. Orientation does not matter
            %   Outputs:    - weight_list:  a list of weights of size [1, N_x*N_y*N_z]
            %               - coor_list:    coor_list is a list of coordinates of size [3, N_x*N_y*N_z]
            % Note: a numerically more robust implementation would probably first do the inner sum, than weigh it with the weights of the middle sum, etc.
            %       our implementation instead pulls the weights together at once.

            %Define lengths and placeholds
            Nx = length(w_x);
            Ny = length(w_y);
            Nz = length(w_z);
            coor_list = zeros(3,Nx*Ny*Nz);
            weight_list = zeros(1,Nx*Ny*Nz);
            
            %Loop through coordinates
            count = 1;
            for i=1:Nx
                for j=1:Ny
                    for k=1:Nz
                        coor_list(:,count) = [coor_x(i); coor_y(j); coor_z(k)];
                        weight_list(1,count) = w_x(i)*w_y(j)*w_z(k);
                        count = count+1;
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform the spatial integration and compute the covariance matrices %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fnc_comp_covariance(obj, interferer_location, target_location, receiver_location) 
            % This function computes the spatial interference and target covariance matrix
            % The correlations between the receiver locations are computed and stored in two covariance matrices (see outputs)
            %   Input:      interferer location (N_interferer x 3), target location (N_target x 3), receiver_location (N_receiver x 3)
            %   Output:     obj.R_Int, obj.R_tar. Both are N_receiver x N_receiver x N_freq_bins

            %We compute the covariance between each combination of microphones. x represents these microphone locations []
            x = receiver_location;
     
            %%%%%%%%%%%%%%%%%%%%%%%
            % Compute covariances %
            %%%%%%%%%%%%%%%%%%%%%%%
            % (1) Interferer
            R_L_INT = 0;
            for i=1:size(interferer_location,1)             % Each loudspeaker is an interferer
                mu = interferer_location(i,:);              % For integration: set mean to loudspeaker location
                R_h = obj.fnc_comp_Rregion(x, mu, obj.IntWinSigma2);    % compute the covariance over the regions using clenshaw curtis quadrature
                R_L_INT = R_L_INT + R_h;                    % We sum the contribution of each region to get the total region
            end
            obj.R_Int = R_L_INT;

            % (2) Target 
            R_L_TAR = 0;
            for i=1:size(target_location,1)
                mu = target_location(i,:);                  % Each person is a target; So far I only used one person so it might not work with more than one.
                R_h = obj.fnc_comp_Rregion(x, mu, obj.TarWinSigma2); % Compute the covariance over the regions.                  
                R_L_TAR = R_L_TAR+R_h;                      % We sum the contribution of each region to get the total region
            end 
            obj.R_Tar = R_L_TAR;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % These are functions which relate to the microphone weights %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fnc_comp_weights(obj, interferer_location, target_location, receiver_location, nyquist_flag)
            %Descr:     This function computes the weights of the microphone spotformer. 
            %
            %Inputs:    R_TAR in Nr x Nr x Nk:     the covariance matrix of the target locations
            %           R_INT in Nr x Nr x Nk:      the covariance matrix of the interferer locations
            %           settings                        the setting object; we use settings.flag_full_axis and settings.nyquist
            %Outputs:   w in Nr x Nk:                   the beamforming weights. 
          
            if nargin == 4
                nyquist_flag = false;
            end
            
            % Get the matrices Rn, Rnum, Riso, R_Int and R_Tar
            Nr = size(receiver_location, 1); %number of receivers
            fnc_comp_Rn(obj, Nr)
            fnc_comp_Risotropic(obj, receiver_location)
            fnc_comp_covariance(obj, interferer_location, target_location, receiver_location) 

            % Get the matrices R_INT = R_Int + R_n + R_num + R_iso; and R_TAR = R_Tar + R_num (numerical inaccuracy is present for both)
            R_INT = obj.R_Int + obj.Rnum + obj.Rn + obj.rebRatio*pagenorm(obj.R_Int).*obj.Riso;
            R_TAR = obj.R_Tar + obj.Rnum;

            disp("Spatially integrating over the regions. This might take a while...")
            for k=1:obj.N_k
                R_INTk = R_INT(:,:,k);    %compute total interferer covariance by summing over the individual ones
                R_TARk = R_TAR(:,:,k);    %idem      
                [V,D] = eig(R_INTk, R_TARk);        %perform generalised eigenvalue decomposition. 
                d = diag(D);
                [min_d,i] = min(d);
                
                %We expect that the imaginary part of d is virtually zero. Similarly, its eigenvalue should be larger than 0 (postive definite)
                if(sum(abs(imag(d)))~=0)
                    disp("Eigenvalues of the " + num2str(k) + "th bin are partially imaginary." ...
                    + "Value (sum of absolute imagninary parts)=" + num2str(sum(abs(imag(d)))))
                end
                if min_d < 0
                    disp("Numerical inaccuracy (I think): the smallest eigenvalue of bin "...
                        + num2str(k) + " is below zero. Value: " +num2str(min_d) );
                end
        
                %Take eigenvector of length 1 corresponding to smallest eigenvalue.
                v = V(:,i);
                w(:,k) = v/norm(v);
                lambda(k) = min_d;
            end
        
            %Artifically set Nyquist bin to zero in case nyquist_flag = true
            if nyquist_flag 
                if obj.flag_full_axis
                    w(:,obj.N_k/2+1) = 0;
                else
                    w(:,end) = 0;
                end
            end
        
            % Assign weights to object 
            obj.weights_mic = w;
        end

        function outputFrame = fnc_comp_output_frame(obj, inputFrame)
            Nr = size(obj.R_Int,1);
            if size(inputFrame) ~= [obj.N_t, Nr]
                disp("inputFrame has incorrect size")
            end

            % FFT of single frame; windowed and zeropadded
            inputFrameHat = fft(obj.w_analysis.*inputFrame, obj.N_t + obj.N_pad);  

            % Multiply with spotformerweights and take ifft
            if obj.flag_full_axis      
                outputFrameHat = obj.weights_mic'.*inputFrameHat;                                               
                outputFrame = ifft(outputFrameHat);
            else
                outputFrameHat = obj.weights_mic'.*inputFrameHat(1:obj.N_t + 1,:);                                    
                outputFrame = ifft([outputFrameHat; zeros(obj.N_pad - 1, Nr)], 'symmetric');
            end
       
            % Sum over the microphones
            outputFrame = obj.w_synthesis.*sum(outputFrame(1:obj.N_t,:),2);   
        end

        function audioOut = comp_output(obj, audioMixture)
            % This function takes as input the audio mixture and outputs the audio as obtained through the microphone spotformer
            %   The advantage of this function over comp_output_headless is that there is no need to recompute the weights for a different audio_input.
            %   However, headless mode does not require first calling the fnc_comp_weights method.
            %
            % Input:     audioMixture   - the LENGTH x Number Recievers audio as obtained by the microphones   
            % Outputs:   audioOut       - the LENGTH x 1 audio signal as obtained through the microphone spotformer. 

            l = 0;              % Frame counter   
            stop_flag = 0;      % Flag indicating when we are done
            audioOut = zeros(size(audioMixture,1)+obj.N_pad, 1); % Set variable containing the output audio     
            while stop_flag == 0
                try %When we are out of frames, an error will be thrown and the stop_flag will be set to 1.
                    inputFrame = audioMixture(l*obj.N_hop+1:l*obj.N_hop+obj.N_t,:); 
                    outputFrame = obj.fnc_comp_output_frame(inputFrame);
                    audioOut(l*obj.N_hop+1:l*obj.N_hop+obj.N_t,:) = audioOut(l*obj.N_hop+1:l*obj.N_hop+obj.N_t,:) + outputFrame;
                catch ME    
                    if strncmp("Index", ME.message, 5)
                        disp("Crashed at frame " + num2str(l) + " v.d. approx " + num2str(floor(size(audioMixture,1)/obj.N_hop - 1)))
                        stop_flag = 1;
                    else    
                        rethrow(ME)
                    end
                end
                l = l+1;
            end
        end


        function audioOut = comp_output_headless(obj, audioMixture, interferer_location, target_location, receiver_location, nyquist_flag)
            % This function takes as input the audio mixture and outputs the audio as obtained through the microphone spotformer
            %   The advantage of this function over comp_output is that you dont need to think about precomputing the weighrs.
            %   The disadvantage is that, if you use the spotformer multiple times in succesion, it is more effecient to compute the weights only once
            %
            % Input:    audioMixture            - the LENGTH x Number Recievers audio as obtained by the microphones   
            %           interferer_location     - [m], (N_interferer x 3), the location of the interferer sound signal
            %           target location         - [m], (N_target x 3), the location of the target sound signal 
            %           receiver_location       - [m], (N_receiver x 3), the location of the receivers (microphones)
            % Outputs:  audioOut                - the LENGTH x 1 audio signal as obtained through the microphone spotformer. 
            if nargin == 5
                nyquist_flag = false;
            end
            obj.fnc_comp_weights(interferer_location, target_location, receiver_location, nyquist_flag)
            audioOut = obj.comp_output(audioMixture);
        end
    end
end

