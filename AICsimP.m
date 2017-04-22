classdef AICsimP<AICBenchP
    %
    properties
        impulseResponse = [];
    end
    methods
       function obj = AICsimP(M,N)
           obj = obj@AICBenchP(M,N);
       end
       function calibrate(obj)
           N = obj.InputRecordLength;
           M = obj.OutputRecordLength;
                    
                T=5e-9; %5 ns
                fs=1e12; % 1 THz
                fp = 200e6; %200 MHz
                f=2e6; %2 MHz
                campioni=T*fs; %Segnale
                dim_prbs=length(obj.prbs)
                livelli=2^(16);
                estr=1;
                
           load Filterresponsenew_19Mhz_800mV
                H=zeros(1,dim_prbs);
                H(1:257)=H0(1:257);
                H(258:dim_prbs-255)=0;
                H(dim_prbs-254:dim_prbs)=conj(fliplr(H0(2:256)));
                obj.impulseResponse = H;
           
            % Direct evaluation of the Phi matrix, fixing the PRBS
            msg = sprintf('Calibration in progress (0/%d)',obj.InputRecordLength);
            HH = waitbar(0,msg);
             for U=0:obj.InputRecordLength/2,
                X(:,U+1)=cos(2*pi*U*(0:obj.InputRecordLength-1)'/obj.InputRecordLength);
             end
              for U=1:(obj.InputRecordLength/2-1),
                X(:,(obj.InputRecordLength/2)+1+U)=sin(2*pi*U*(0:obj.InputRecordLength-1)'/obj.InputRecordLength);
              end
            for J=0:obj.InputRecordLength/2,
                X2=cos(2*pi*J*(0:dim_prbs-1)'/dim_prbs);
                Y(:,J+1)=generateAndAcquire(obj,X2);
%                 mout = X2.*obj.prbs;
%                 MOUT=fft(mout);
%                 SIGF=real(ifft((MOUT.').*H));
%                 Y(:,J+1)=SIGF(2500:5000*obj.undersamplingRate:end); 
%                 partition = [-estr:(1/livelli):estr];
%                 codebook = [(-estr-1/livelli):(1/livelli):estr];
%                 [index,]=quantiz(SIGF2,partition,codebook);
                msg = sprintf('Calibration in progress (%d/%d)',J,obj.InputRecordLength);
                waitbar(J/obj.InputRecordLength,HH,msg);
%                 size(mout)
%                 size(MOUT)
%                 size(H)
%                 size(ifft((MOUT)))
%                 LUGHEZZASIGF=size(SIGF)
            end
            for J=1:(obj.InputRecordLength/2-1),
                X2=sin(2*pi*J*(0:dim_prbs-1)'/dim_prbs);
                Y(:,(obj.InputRecordLength/2)+1+J)=generateAndAcquire(obj,X2);
%                 mout = X2.*obj.prbs;
%                 MOUT=fft(mout);
%                 SIGF=real(ifft((MOUT.').*H));
%                   partition = [-estr:(1/livelli):estr];
%                   codebook = [(-estr-1/livelli):(1/livelli):estr];
%                 Y(:,(obj.InputRecordLength/2)+1+J)=SIGF(2500:5000*obj.undersamplingRate:end);
%                 [index,]=quantiz(SIGF2,partition,codebook);
                msg = sprintf('Calibration in progress (%d/%d)',J+obj.InputRecordLength/2+1,obj.InputRecordLength);
                waitbar((obj.InputRecordLength/2+J+1)/obj.InputRecordLength,HH,msg);
            end
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
           sigout2 = SIGF(2500:obj.undersamplingRate*5000:end);
           livelli=2^(16);
           estr=2;
           x=linspace(-estr,estr,livelli);
           sigout=interp1(x,x,sigout2,'nearest');
%            sigout=floor((sigout2.*livelli)./estr)./livelli.*estr;
%            partition = [-estr:(1/livelli):estr];
%            codebook = [(-estr-1/livelli):(1/livelli):estr];
%            [index,sigout]=quantiz(sigout2,partition,codebook);
%            Yq=round(sigout2/estr*livelli/2-1)/(livelli/2-1)*estr;
%            set=find(abs(Yq)>estr);
%            Yq(set)=estr*sign(Yq(set)); 
%            sigout=Yq;
       end
    end
end