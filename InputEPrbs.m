% testAICsim
% for iter=1:2,
clc
clear all
% clearvars -except prbsa fasecasuale
% clearvars -except iter
% close all
N=512;
M = 128;
eM=4;
T=5e-9; %5 ns
fs=1e12; % 1 THz
fp = 200e6; %200 MHz
f=200; %2 MHz
campioni=T*fs; %Segnale
dim_prbs=N*campioni;

bench = AICsimPC(N,M); %crea ogetto AICsim->Funzioni: Calibrate; generateAndAcquire


% prbsa=randi(2,N+2,1)-1; %dirac pulse +2

% i=1;
% for jitt=0
% for W=100:10:100
% for W=0:100:300;
W=100;
figure;
jitt=4;
        for randomj=1:1;
            prbsa=randi(2,N+2,1)-1; %dirac pulse +2
            
        bench.init(prbsa, jitt ,W); %crea la prbs all'ogetto
        %0 NoJitt, 1 DCD, 2 ISI, 3 PJ, 4 RJ
        %W=0:10:300; per 1 e 3
        %W2=0:10:100; per 2 e 4
        Normaliz=sum(abs(bench.prbs))/dim_prbs;
        disp('Initialization done');
        
        disp('Starting calibration...');
        bench.calibrate(); %dichiara filtro e crea measurementMatrix all'ogetto

        % INGRESSO
        freqVec = [125]'; 
        ampVec = [0]'; 
%         fasecasuale=rand(length(freqVec),1);
%%      
%         figure;
        i=1;
        hold on;
%         figure;
        fasemisurate=[ .5 ];
%         fasemisurate=[0 .1 .5 .9 1];
        for f=1:length(fasemisurate);
            phi = 2*pi*fasemisurate(f);
            
            X = zeros(dim_prbs,1);
            X(freqVec+1) = (10.^(ampVec/20)).*exp(1i*phi)./Normaliz;
            X(dim_prbs/2+2:dim_prbs) = conj(flipud(X(2:dim_prbs/2))); 
            x = ifft(X)*dim_prbs; %reale
                X2=zeros(N,1);
                X2(freqVec+1)=X(freqVec+1);
                X2(N/2+2:N) = conj(flipud(X2(2:N/2))); 
                x2=ifft(X2)*N;
            y(i,:) = bench.generateAndAcquire(x); %moltiplica prbs e filtra
    %         figure;
    %         plot(y(i,:));
    %         hold on;
    %         y2(:,i)=bench.measurementMatrix*x2;
    %         plot(bench.measurementMatrix*x2,'r');
    %         hold off
    %         diff(:,i)=y(i,:)'-y2(:,i);
    %         figure; plot(y(i,:)'-bench.measurementMatrix*x2)

            epsilon=[1e-8 1.5573827825e-8 2e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10];
%             epsilon=[2.5e-16 1e-15 1e-14 1e-13 1e-12 1e-11  1e-10  1e-9 1e-8 2e-8 1e-7 1e-6 ];
%             epsilon=[3e-16 1e-15 1e-14 1e-13 1e-12 1e-11  1e-10  1e-9 1e-8 2e-8 1e-7 1e-6 ];
            for k=1:length(epsilon)
                lun=length(epsilon)*(i-1);
                [xest,xp]=bench.reconstruct(y(i,:)','Dantzig',epsilon(k)); % crea Psi e moltiplica per measurementMatrix poi l1magic
                Xest(:,k+lun) = fft(xest)/N;
%                     figure(2*iter);
%                     plot(0:N/2,db(abs(Xest(1:N/2+1,k+lun))));
%                     hold on;
%                     plot(freqVec,ampVec,'or');
%                     hold off
                Xest(1,k+lun)=0;
                maxXest=find(abs(Xest(:,k+lun))==max(abs(Xest(:,k+lun))));
                Xest2=Xest(:,k+lun);
                Xest2(maxXest(1),1)=0;
                Xest2(maxXest(2),1)=0;
                maxXest2=find(abs(Xest2(:,1))==max(abs(Xest2(:,1))),2);

                SFDR(i,k)=20*log10((abs(Xest(maxXest(1),k+lun)))/(abs(Xest2(maxXest2(1),1))));
                SINAD(i,k)=20*log10(sqrt(((N-3)/N)*((2*((abs(Xest(maxXest(1),k+lun))).^2))/(sum(abs(Xest2(:,1)).^2)))));
            end
%             figure(2*iter+1);
            hold on
            semilogx(epsilon(:),SINAD(i,:),'-o');
            set(gca,'xscale','log')
            xlabel('Epsilon')
            ylabel('SINAD')
            grid on
            hold off
%             legendInfo{i} = ['Fase = ' num2str(fasemisurate(f))];
            i=i+1;
%          end
        end
%     legend(legendInfo);
        end
    title(['Random Jitter per W=' num2str(W)]);
% end
% end

