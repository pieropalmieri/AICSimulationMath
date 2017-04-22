classdef AICBenchP<handle
    properties
        prbsRate;
        prbs = [];
        
        InputRecordLength;
        OutputRecordLength;
        undersamplingRate;
        sampleRate;
        measurementMatrix = [];
    end
    
    methods
        function obj = AICBenchP(N,M)
            obj.InputRecordLength = N;
            obj.OutputRecordLength = M;
            obj.undersamplingRate = N/M;
            obj.sampleRate = obj.prbsRate/obj.undersamplingRate;
            
            obj.prbs = zeros(N,1);
        end
        
        function init(obj, prbsa, jitt, W)

            if ~exist('jitt','var')
                jitt=0;
            end
            if ~exist('W','var')
                W=0;
            end
            
                N=obj.InputRecordLength;
                eM=obj.undersamplingRate;
                M=obj.OutputRecordLength;

%                 prbsa=randi(2,N+2,1)-1; %dirac pulse +2
                
                T=5e-9; %5 ns
                fs=1e12; % 1 THz
                campioni=T*fs; %Segnale
                transizione=625;

                load -ascii ISITrace00000_old.dat;
                R0=reshape(ISITrace00000_old,250,length(ISITrace00000_old)/250);
                R=zeros(size(R0,2),campioni+transizione*2); %Campioni + ?

                for index=1:size(R0,2)
                    R(index,:)=interp1(-625e-12:25e-12:25e-12*199+625e-12,R0(:,index),-625e-12:1e-12:4999e-12+625e-12);
                end % da -625ps con passo 25ps a passo 1ps

                %transizioni 010, 011, 110
                R(3,1+155:3125)=R(3,1:3125-155); %shift di 155 avanti
                R(3,3126:6250-55)=R(3,3126+55:6250); %shift di 55 valori indietro
                R(4,1+148:6250)=R(4,1:6250-148);
                R(7,1:6250-76)=R(7,1+76:6250);
                %transizioni 101, 001, 100
                R(6,1:3125-89)=R(6,89+1:3125);
                R(6,3125-88:3125+182)=R(6,3125);
                R(6,3126+182:6250)=R(6,3126:6250-182);
                R(2,1+160:6250)=R(2,1:6250-160);
                R(5,1:6250-80)=R(5,80+1:6250);
                
                
                if(jitt==1)
                    DCD=[0;-1;1;0;0;-1;1;0]*W;
                    prbs_DCD=zeros(N*campioni,1);
                    r_DCD=zeros(length(W),1);
                    for i=1:length(W)
                        conta=0;
                        for n=2:N+1
                            index=binaryVectorToDecimal(prbsa(n-1:n+1)')+1;
                            prbs_DCD(1+conta:campioni/2+conta)=R(index,transizione+1-DCD(index,i):transizione+campioni/2-DCD(index,i));
                            prbs_DCD(campioni/2+1+conta:campioni+conta)=R(index,transizione+campioni/2+1-DCD(index,i):transizione+campioni-DCD(index,i));
                            conta=conta+campioni;
                        end
                        ampiezza=(49.10-(-47.79))/2/1000;
                        prbs_DCD=1/ampiezza*prbs_DCD;
                        obj.prbs=prbs_DCD;
                   end
                    disp('Hai scelto DCD')
                elseif(jitt==2)
                        ISI=[0;1;-1;0;0;-1;1;0]*W;
                        prbs_ISI=zeros(N*campioni,1);
                        r_ISI=zeros(length(W),1);
                        for i=1:length(W)
                            conta=0;
                            for n=2:N+1
                                index=binaryVectorToDecimal(prbsa(n-1:n+1)')+1;
                                prbs_ISI(1+conta:campioni/2+conta)=R(index,transizione+1-ISI(index,i):transizione+campioni/2-ISI(index,i));
                                prbs_ISI(campioni/2+1+conta:campioni+conta)=R(index,transizione+campioni/2+1-ISI(index,i):transizione+campioni-ISI(index,i));
                                conta=conta+campioni;
                            end
                            ampiezza=(49.10-(-47.79))/2/1000;
                            prbs_ISI=1/ampiezza*prbs_ISI;
                             obj.prbs=prbs_ISI;
                        end
                        disp('Hai scelto ISI')
                  elseif(jitt==3)
                        PJ=W;
                        prbs_PJ=zeros(N*campioni,1);
                        r_PJ=zeros(length(W),1);
                        for i=1:length(W)
                            conta=0;
                            for n=2:N+1
                                index=binaryVectorToDecimal(prbsa(n-1:n+1)')+1;
                                prbs_PJ(1+conta:campioni/2+conta)=R(index,transizione+1-PJ(i):transizione+campioni/2-PJ(i));
                                prbs_PJ(campioni/2+1+conta:campioni+conta)=R(index,transizione+campioni/2+1-PJ(i):transizione+campioni-PJ(i));
                                conta=conta+campioni;
                                PJ(i)=-PJ(i);
                            end
                            ampiezza=(49.10-(-47.79))/2/1000;
                            prbs_PJ=1/ampiezza*prbs_PJ;
                             obj.prbs=prbs_PJ;
                        end
                        disp('Hai scelto PJ')
                 elseif(jitt==4)
                  sigma=W'*1e-12;
%                     iter=20;
%                     iter=1;
                    prbs_rj=zeros(N*campioni,1);
%                     r_rj=zeros(length(sigma),iter);
%                     r_rj=zeros(length(sigma));
%                     for j=1:iter
                        gauss=randn(N+1,1);
%                         for k=1:length(sigma)
%                             tie=sigma(k)*gauss;
                            tie=sigma*gauss;
                            campioni_rj=mod(round(tie*fs),transizione);
%                             max(abs(campioni_rj))
                            conta=0;
                            for n=2:N+1
                                index=binaryVectorToDecimal(prbsa(n-1:n+1)')+1;
                                prbs_rj(1+conta:campioni/2+conta)=R(index,transizione+1-campioni_rj(n-1):transizione+campioni/2-campioni_rj(n-1));
                                prbs_rj(campioni/2+1+conta:campioni+conta)=R(index,transizione+campioni/2+1-campioni_rj(n):transizione+campioni-campioni_rj(n));
                                conta=conta+campioni;
                            end
                            ampiezza=(49.10-(-47.79))/2/1000;
                            prbs_rj=1/ampiezza*prbs_rj;
                            obj.prbs=prbs_rj;
%                         end
%                     end
                  disp('Hai scelto RJ')
                else
                prbs_nojitter=zeros(N*campioni,1); %dim_prbs
                conta=0;
                for n=2:N+1
                    index=binaryVectorToDecimal(prbsa(n-1:n+1)')+1; %3 bit fino a 8
                    prbs_nojitter(1+conta:campioni/2+conta)=R(index,transizione+1:transizione+campioni/2);
                    prbs_nojitter(campioni/2+1+conta:campioni+conta)=R(index,transizione+campioni/2+1:transizione+campioni);
                    conta=conta+campioni;
                end 
                
                ampiezza=(49.10-(-47.79))/2/1000;
                prbs_nojitter=1/ampiezza*prbs_nojitter;
                
%                 prbs_nojitterD = (randi(2,N,1)-1)*2-1;
%                 obj.prbs=zeros(campioni*N,1);
%                 obj.prbs(1:campioni:campioni*N,1)=prbs_nojitterD(1:end,1);
                
                obj.prbs=prbs_nojitter;
%                 obj.prbs=prbs_nojitterD;
                disp('Hai scelto NO jitter')
                end
               fprintf('Il W scelto è: %d \n',W)
        end
        
        function [sigest,xp]=reconstruct(obj,aicout,algorithm,par)
            
            N = obj.InputRecordLength;
            
            Psi = (cos(2*pi*(0:N-1)'*(0:N-1)/N)+sin(2*pi*(0:N-1)'*(0:N-1)/N))/sqrt(N);
            %DHT da vedere la trasformata
            A = obj.measurementMatrix * Psi;
            
            if strcmp(algorithm,'Dantzig'),
                % L1-magic Dantzig
                w = A*A' \ aicout;
                x0 = A'*w;

                % Dantzig selection

                xp = l1dantzig_pd(x0, A, [], aicout, par, eps);
            elseif strcmp(algorithm,'CoSaMP'),
                xp = CoSaMP( A, aicout, par, [],[]);
            end
            sigest = Psi * xp;
        end
    end
    
    methods (Abstract)
        calibrate(obj, prbsa, jitt, W)
        sigout=generateAndAcquire(obj,sigin)
    end  
        
end