classdef AICsimPC<AICBenchP
    
    properties
        impulseResponse = [];
    end
    methods
       function obj = AICsimPC(M,N)
           obj = obj@AICBenchP(M,N);
       end
       function calibrate(obj)
           N = obj.InputRecordLength;
           M = obj.OutputRecordLength;
           
           %obj.measurementMatrix = randn(M,N);
           
%            L = 32;
%            h = fir1(L,1/5);
%            obj.impulseResponse = h;
                fp = 200e6; %200 MHz
                T=5e-9; %5 ns
                fs=1e12; % 1 THz
                campioni=T*fs; %Segnale
                dim_prbs=N*campioni
%                 dim_prbs=length(obj.prbs)
           
           load Filterresponsenew_19Mhz_800mV
                 H=zeros(1,dim_prbs);
                H(1:257)=H0(1:257);
                H(258:dim_prbs-255)=0;
                H(dim_prbs-254:dim_prbs)=conj(fliplr(H0(2:256)));
                obj.impulseResponse = H;
                
                K=fp/fs;
                L0=((fp+fs)/(2*fp/N))-1;
                L=2*L0+1
                s=obj.prbs;
                length(obj.prbs);
                for p=1:(L),
                    for q=1:(dim_prbs),
                    F(q,1)=exp((-1i*2*pi/N)*p*q);
%                     (q/dim_prbs)*100
                    end
                    molt(p)=s*F(q,1);
                    size(molt)
                end
%                 size(G)
%                 for p=1:L,
%                     F(1,p)=G(:,p);
%                 end
                for p=1:L,
                    if p==0,
                        d(p)=1/N;
                    else
                        d(p)=(1-exp((-1i*2*pi/N)*p))/(2*pi*1i*p);
                    end
                end
                D=diag(d);
                Phi=s*F*D;
                obj.measurementMatrix=Phi;

%            
%             % Direct evaluation of the Phi matrix, fixing the PRBS
%             msg = sprintf('Calibration in progress (0/%d)',dim_prbs);
%             HH = waitbar(0,msg);
%             for J=0:dim_prbs/2,
%                 X(:,J+1)=cos(2*pi*J*(0:dim_prbs-1)'/dim_prbs);
%                 mout = X(:,J+1).*obj.prbs;
%                 MOUT=fft(mout);
%                 SIGF=real(ifft((MOUT.').*H));
%                 Y(:,J+1)=SIGF(1:obj.undersamplingRate:length(SIGF));
%                 
% %                 sigf = filter(obj.impulseResponse,1,mout);
% %                 Y(:,J+1) = sigf(1:obj.undersamplingRate:length(sigf));
% %                 
%                 msg = sprintf('Calibration in progress (%d/%d)',J,dim_prbs);
%                 waitbar(J/dim_prbs,HH,msg);
%             end
%             for J=1:dim_prbs/2,
%                 X(:,(dim_prbs+1)/2+J)=sin(2*pi*J*(0:dim_prbs-1)'/dim_prbs);
%                 
%                 mout = X(:,(dim_prbs+1)/2+J).*obj.prbs;
%                 MOUT=fft(mout);
%                 SIGF=real(ifft((MOUT.').*H));
%                 Y(:,(dim_prbs+1)/2+J)=SIGF(1:obj.undersamplingRate:length(SIGF));
%                 
% %                 sigf = filter(obj.impulseResponse,1,mout);
% %                 Y(:,N/2+1+J) = sigf(1:obj.undersamplingRate:length(sigf));
% %                 
%                 msg = sprintf('Calibration in progress (%d/%d)',J+dim_prbs/2+1,dim_prbs);
%                 waitbar((J+dim_prbs/2+1)/dim_prbs,HH,msg);
%             end
            close(HH);
            size(Y)
            size(X)
            Phi=Y*inv(X);
            obj.measurementMatrix = Phi;

% % %            H = toeplitz([h'; zeros(N-L,1)],[h(1) zeros(1,N-1)]);
% % %            D = diag(obj.prbs);
% % %            
% % %            PhiEx = H*D;
% % %            obj.measurementMatrix = PhiEx(1:obj.undersamplingRate:N,:);
       end
       
       function sigout=generateAndAcquire(obj,sigin)
           %sigout = obj.measurementMatrix*sigin;
           mout = sigin.*obj.prbs;
           MOUT=fft(mout);
           SIGF=real(ifft((MOUT.').* obj.impulseResponse ));
           sigout = SIGF(2500:obj.undersamplingRate*5000:end);
%            H0=obj.impulseResponse;
%                 H=zeros(1,dim_prbs);
%                 H(1:257)=H0(1:257);
% %                 H(258:dim_prbs-255)=0; % ?
%                 H(dim_prbs-254:dim_prbs)=conj(fliplr(H0(2:256)));
%                 H(258:512)=conj(fliplr(H0(2:256)));
%            SIGF=real(ifft((MOUT.').* obj.impulseResponse));
% %            sigf = filter(obj.impulseResponse,1,mout);
%            sigout = SIGF(1:obj.undersamplingRate:length(SIGF));
       end
    end
end