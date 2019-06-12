clear; clc; close all;

%Get points
path = '.\54_HOA2-HA-E21211-HA-2s\20190319_083641';     %GN=24.0dB

saveFlag = 0;

nomGain = 24;   %Nominal Gain
stg = str2num(path(strfind(path,'s\')-1));


% mask
data = dlmread([path '\mask\mask.txt'],'	',2,0);
    targetGain       = data(:,1);
    pin              = data(:,2);
    targetOutput     = data(:,3);
    laser3Pwr        = data(:,11);
    srcCh            = data(:,14);
    ampCh            = data(:,15);


%Spectrum
files = dir([path '\spectrum\']);
count = 0;
for j = 3:length(files)
    if any(findstr(files(j).name,'Src'))  
        continue   
    end
    count = count+1;  %pega o conteudo de Amp.  
    id = findstr(files(j).name,'_');
    TargetGain(count) = str2num(files(j).name(id(1)+1:id(2)-1));
    Pin(count)        = str2num(files(j).name(id(3)+1:id(4)-1));
end
TargetPout = Pin + TargetGain;

%Sort
aux = sortrows([TargetGain' Pin' TargetPout'],[1 2]);   % Ordena como? imprimir!
TargetGain  = aux(:,1)';
Pin         = aux(:,2)';
TargetPout  = aux(:,3)';
    

planck      = 6.626068e-34;
light       = 299792458;
res         = 0.097 * 1e-9; %0.1 * 1e-9

%IL
if stg == 1
    %STG 1
    Fiber = 19.96;
elseif stg == 2
    %STG 2
    Fiber = 0;
end


NoiseConv = 0;  %10*log10(0.3/(res*1e9));  %Converter o nível de ruído para a banda do sinal (37.5GHz)


%AllLamS
spectrumSrc = dlmread([[path '\spectrum\'] sprintf('Gain_%05.2f_Pin_%05.2f_Src.txt',max(TargetGain)-2,max(TargetPout)-max(TargetGain))],'	',2,0);
Dres = (max(spectrumSrc(:,1))-min(spectrumSrc(:,1)))/(length(spectrumSrc)-1);
[pmaxSrc,loc] = findpeaks(spectrumSrc(:,2),'MinPeakHeight',-60,'MinPeakProminence',20,'MinPeakDistance',ceil(3*res/Dres));
AllLamS = spectrumSrc(loc,1)'*1e9;
AllLamS = round(AllLamS*10)/10;


%IL
Src2UUT = 0;
Src2PM = [];
UUT2PM = [];
ramanGain = [];
ErrorPoint= [];
PoutPM = [];
PoutPD = [];

flag = 1;
for j = 1:length(TargetGain)
    pData = dlmread([[path '\result\'] sprintf('Gain_%05.2f_Pin_%05.2f_PM.txt',TargetGain(j),Pin(j))],' ',1,0);
    PoutPM = [PoutPM pData(2)];
    if stg == 1
        %STG 1
        PoutPD = [PoutPD pData(4)];
    elseif stg == 2
        %STG 2
        PoutPD = [PoutPD pData(6)];
    end
    
    
    if pData(9)>-100
        ramanGain = [ramanGain; pData(3)-pData(9)];
    elseif ~isempty(ramanGain)
        ramanGain = [ramanGain; ramanGain(end)];
    end
    
%     if TargetGain(j)==38 && TargetPout(j)==19
%         keyboard;
%     end
    
    if stg == 1
        %STG 1
        if pData(2)<-100 || pData(9)<-45
            continue
        end
        Src2PM = [Src2PM; (pData(9)+Fiber)-pData(1)];
        UUT2PM = [UUT2PM; pData(4)-pData(2)];
    elseif stg == 2
        %STG 2
        if pData(2)<-100
            continue
        end
        Src2PM = [Src2PM; pData(5)-pData(1)];
        UUT2PM = [UUT2PM; pData(6)-pData(2)];
    end
    
%     if Pin(j)==max(Pin)-10 && TargetGain(j)==max(TargetGain)-5
    if TargetGain(j)==max(TargetGain)-8 && TargetPout(j)>max(TargetPout)-9 && flag>0
        flag = 0;
        spectrumSrc = dlmread([[path '\spectrum\'] sprintf('Gain_%05.2f_Pin_%05.2f_Src.txt',TargetGain(j),Pin(j))],'	',2,0);
        factor  = 0.000215 * length(spectrumSrc)/(0.1);
        power  = 10*log10(trapz(10.^(spectrumSrc(:,2)/10)*1e-3)*1/factor*1e3);
        OSA2PM = pData(1)-power;
        if length(pData)>6
            %(1)PMIn	(2)PMOut	(3)PDIn Stg1	(4)PDOut Stg1	(5)PDIn Stg2	(6)PDOut Stg2	(7)VOA Stg1	(8)VOA Stg2	(9)PDIn OFF
            if stg == 1
                %STG 1
                nomGain = round((TargetGain(j)-(pData(3)-pData(9))+pData(7))*10)/10;    %TargetGain - RamanGain + VOA
            elseif stg == 2
                %STG 2
                nomGain = round((TargetGain(j)+pData(8))*10)/10;
            end
        end
    end

end
Src2PM = mean(Src2PM);
UUT2PM = mean(UUT2PM);

fprintf('Src2PM: %1.2f dB\nUUT2PM: %1.2f dB\nOSA2PM: %1.2f dB\n\n',Src2PM,UUT2PM,OSA2PM);

PoutPM = PoutPM + UUT2PM;


% iteration
chGain      = zeros(length(TargetGain),length(AllLamS));
chNF        = zeros(length(TargetGain),length(AllLamS));
for j = 1:length(TargetGain)
        
        spectrumSrc = dlmread([[path '\spectrum\'] sprintf('Gain_%05.2f_Pin_%05.2f_Src.txt',TargetGain(j),Pin(j))],'	',2,0);
        spectrumAmp = dlmread([[path '\spectrum\'] sprintf('Gain_%05.2f_Pin_%05.2f_Amp.txt',TargetGain(j),Pin(j))],'	',2,0);
        
        spectrumSrc(:,2) = spectrumSrc(:,2) + (Src2PM - Src2UUT + OSA2PM);
        spectrumAmp(:,2) = spectrumAmp(:,2) + (UUT2PM + OSA2PM);
        
        
        if TargetGain(j)==28 && TargetPout(j)==-7
            figure('windowstyle','docked'); hold on; grid on; box on;
            plot(spectrumSrc(:,1),spectrumSrc(:,2),'LineWidth',1)
            plot(spectrumAmp(:,1),spectrumAmp(:,2),'LineWidth',2)
            xlabel('\lambda (nm)')
            ylabel('Optical Power (dBm)')
            title(sprintf('(Target %d dBm; Pout %d dBm)',TargetGain(j),TargetPout(j)))
            set(gca, 'FontSize',20)
        end
        
        
        factor  = 0.000215 * length(spectrumSrc)/(0.1);
        PowerIN(j)  = 10*log10(trapz(10.^(spectrumSrc(:,2)/10)*1e-3)*1/factor*1e3);
        PowerOUT(j) = 10*log10(trapz(10.^(spectrumAmp(:,2)/10)*1e-3)*1/factor*1e3);
        RealTotalGain(j) = PowerOUT(j) - (PowerIN(j) - Fiber);
        
        
        if abs(RealTotalGain(j)-TargetGain(j))>5   %AGCErrorEval
            fprintf('NOK AGC: Gain_%05.2f_Pin_%05.2f.txt\n',TargetGain(j),Pin(j));
            ErrorPoint = [ErrorPoint; TargetGain(j) Pin(j)];
        end

        if min(PowerIN(j),PowerOUT(j))<-50
            fprintf('NOK SPECTRUM: Gain_%05.2f_Pin_%05.2f.txt\n',TargetGain(j),Pin(j));
            ErrorPoint = [ErrorPoint; TargetGain(j) Pin(j)];
            continue
        end
        
        
        Dres = (max(spectrumAmp(:,1))-min(spectrumAmp(:,1)))/(length(spectrumAmp)-1);
        
        %MinPeakHeight: -65(stg1) ou -50(stg2)
        [pmaxSrc,loc] = findpeaks(spectrumSrc(:,2),'MinPeakHeight',-50,'MinPeakProminence',20,'MinPeakDistance',ceil(3*res/Dres));
        pmaxSrc = [spectrumSrc(loc,1)*1e9 pmaxSrc];
                
        sltRange = spectrumAmp(:,1)*1e9>min(pmaxSrc(:,1))-7 & spectrumAmp(:,1)*1e9<max(pmaxSrc(:,1))+7;
        sltSpectrumAmp = spectrumAmp(sltRange,1);
        sltSpectrumAmp = [sltSpectrumAmp spectrumAmp(sltRange,2)];
        MinPeak = max(sltSpectrumAmp(:,2))-20;
        [pmaxAmp,locMax] = findpeaks(sltSpectrumAmp(:,2),'MinPeakHeight',MinPeak,'MinPeakProminence',12,'MinPeakDistance',ceil(3*res/Dres));
        pmaxAmp = [sltSpectrumAmp(locMax,1)*1e9 pmaxAmp];
        
        
        if length(pmaxAmp)==0
            fprintf('NOK PEAK: Gain_%05.2f_Pin_%05.2f.txt\n',TargetGain(j),Pin(j));
            ErrorPoint = [ErrorPoint; TargetGain(j) Pin(j)];
            continue
        end
        
        if any(ampCh(targetGain==TargetGain(j) & pin==Pin(j))) && length(pmaxAmp)<ampCh(targetGain==TargetGain(j) & pin==Pin(j))
            fprintf('NOK OSNR: Gain_%05.2f_Pin_%05.2f.txt\n',TargetGain(j),Pin(j));
            ErrorPoint = [ErrorPoint; TargetGain(j) Pin(j)];
        end
        
        
        loc = [];
        for i = 1:length(locMax)-1
            id = find(sltSpectrumAmp(:,2)==min(sltSpectrumAmp(locMax(i):locMax(i+1),2)));
            id(id<locMax(i) | id>locMax(i+1)) = [];
            loc = [loc; id(1)];
        end
        %extremo esquerdo
        id = find(sltSpectrumAmp(:,2)==min(sltSpectrumAmp(locMax(1)-ceil(5*res/Dres):locMax(1),2)));
        id(id<(locMax(1)-ceil(5*res/Dres)) | id>locMax(1)) = [];
        loc = [loc; id(1)];
        %extremo direito
        id = find(sltSpectrumAmp(:,2)==min(sltSpectrumAmp(locMax(end):locMax(end)+ceil(5*res/Dres),2)));
        id(id<locMax(end) | id>(locMax(end)+ceil(5*res/Dres))) = [];
        loc = [loc; id(1)];
        
        pminAmp = [sltSpectrumAmp(loc,1)*1e9 sltSpectrumAmp(loc,2)];
        
        
        sltSpectrumSrc = spectrumSrc(sltRange,1);
        sltSpectrumSrc = [sltSpectrumSrc spectrumSrc(sltRange,2)];
        pmaxSrc = [sltSpectrumSrc(locMax,1)*1e9 sltSpectrumSrc(locMax,2)];
        pminSrc = [sltSpectrumSrc(loc,1)*1e9 sltSpectrumSrc(loc,2)];
        
        
        InFlatnessEval(j) = max(pmaxSrc(:,2))-min(pmaxSrc(:,2));
        if InFlatnessEval(j) > 1.0
            fprintf('NOK EQ: Gain_%05.2f_Pin_%05.2f.txt\n',TargetGain(j),Pin(j));
            ErrorPoint = [ErrorPoint; TargetGain(j) Pin(j)];
        end

        
        
        lamS        = pmaxAmp(:,1);
        Freq        = light./lamS * 1e9;
        DeltaFreq   = (Freq.^2).*res./light;
        
        
        PaseAuxIN = ([10.^((pminSrc(:,2)+NoiseConv)/10);0] + [0;10.^((pminSrc(:,2)+NoiseConv)/10)])./2;
        PaseAuxOUT= ([10.^((pminAmp(:,2)+NoiseConv)/10);0] + [0;10.^((pminAmp(:,2)+NoiseConv)/10)])./2;
        PaseIN  = [lamS PaseAuxIN(2:end-1)*1e-3];
        PaseOUT = [lamS PaseAuxOUT(2:end-1)*1e-3];
        
        RealChannelGain  = [lamS 10*log10((10.^(pmaxAmp(:,2)/10)*1e-3-PaseOUT(:,2))*1e3) - 10*log10((10.^(pmaxSrc(:,2)/10)*1e-3-PaseIN(:,2))*1e3)];
        
        NF_linear = (PaseOUT(:,2) - PaseIN(:,2).*10.^(RealChannelGain(:,2)/10))./(planck*Freq.*DeltaFreq.*10.^(RealChannelGain(:,2)/10));
        NFsn_linear = NF_linear + 1./10.^(RealChannelGain(:,2)/10);
        NF_linear(NF_linear<0) = NaN;
        NFsn_linear(NFsn_linear<0) = NaN;

        NFEval(j) = 10*log10(max(NF_linear));
        NFsnEval(j) = 10*log10(max(NFsn_linear));
        
        GainFlatnessEval(j) = max(pmaxAmp(:,2)) - min(pmaxAmp(:,2));
        

        id = [];
        for k = 1:length(lamS)
            id = [id find(abs(round(pmaxSrc(k,1)*10)/10-AllLamS)==min(abs(round(pmaxSrc(k,1)*10)/10-AllLamS)))];
        end
        chGain(j,id) = RealChannelGain(:,2) + Fiber;
        chNF(j,id)   = 10*log10(NF_linear) - Fiber;
        pmax    = zeros(1,length(AllLamS)); pmax(id)   = pmaxAmp(:,2);
        chPaseOUT= zeros(1,length(AllLamS));chPaseOUT(id)= 10*log10(PaseOUT(:,2)*1e3);
        chPaseIN =zeros(1,length(AllLamS)); chPaseIN(id) = 10*log10(PaseIN(:,2)*1e3);
        
end

AGCErrorEval = RealTotalGain - TargetGain;
PowerIN = PowerIN - Fiber;
NFEval = NFEval - Fiber;
NFsnEval = NFsnEval - Fiber;



%% Plot
sz = 100;

% Points
%           ________  __PoutMax
%         /|_______ /|__Pout64
%        /_|______ /_|__Pout32
%       /__|_____ /__|__Pout16
%      /___|____ /___|__Pout8
%     /____|___ /____|__Pout4
%    /_____|__ /_____|__Pout2
%    |     |   |     |
%  Pin1  Pin2 Pin3  Pin4

if stg == 1
    % HOA1
    MaxGain = 38;
    IntGain = 35;   %Ganho controlado
    MinGain = 20;
    PoutMax = max(TargetPout);  %20 ou 21 dBm
    Pin1    = -43;
elseif stg == 2
    %EOA1 HG
    MaxGain = 28;
    MinGain = 10;
    PoutMax = 21;
    Pin1    = -35;
end


Pin2    = PoutMax-MaxGain;
Pin4    = PoutMax-MinGain;
Pin3    = Pin1+(Pin4-Pin2);
Pout2   = Pin1+MaxGain;
Pout64  = PoutMax-3;
Pout32  = Pout64-3;
Pout16  = Pout32-3;
Pout8   = Pout16-3;
Pout4   = Pout8-3;
if exist('IntGain')
    PinI1   = Pin1+(MaxGain-IntGain);
    PinI2   = Pin2+(MaxGain-IntGain);
end


figure('windowstyle','docked'); hold on; grid on; box on;
scatter(Pin,TargetPout,sz,InFlatnessEval,'filled'); colorbar;
axis([Pin1-2 PoutMax-MinGain+2 Pin1+MaxGain-2 PoutMax+2])
xlabel('Pin (dBm)')
ylabel('Pout (dBm)')
title('SrcFlatness')
set(gca, 'FontSize',20)


if ~isempty(ErrorPoint)
    sltNOKPoint = ismember([TargetGain' Pin'],ErrorPoint,'rows');
    
    % --------------------------------------------- %
    PoutPM(sltNOKPoint)                        = [];
    PoutPD(sltNOKPoint)                        = [];
    Pin(sltNOKPoint)                           = [];
    PowerOUT(sltNOKPoint)                      = [];
    TargetGain(sltNOKPoint)                    = [];
    TargetPout(sltNOKPoint)                    = [];
    RealTotalGain(sltNOKPoint)                 = [];
    GainFlatnessEval(sltNOKPoint)              = [];
    NFEval(sltNOKPoint)                        = [];
    AGCErrorEval(sltNOKPoint)                  = [];
    % --------------------------------------------- %
end


figure('windowstyle','docked'); hold on; grid on; box on;
scatter(Pin,TargetPout,sz,RealTotalGain,'filled'); colorbar;
axis([Pin1-2 PoutMax-MinGain+2 Pin1+MaxGain-2 PoutMax+2])
xlabel('Pin (dBm)')
ylabel('Pout (dBm)')
title('RealTotalGain')
set(gca, 'FontSize',20)

figure('windowstyle','docked'); hold on; grid on; box on;
scatter(Pin,TargetPout,sz,GainFlatnessEval,'filled'); colorbar;
axis([Pin1-2 PoutMax-MinGain+2 Pin1+MaxGain-2 PoutMax+2])
xlabel('Pin (dBm)')
ylabel('Pout (dBm)')
title('Gain Flatness')
set(gca, 'FontSize',20)

figure('windowstyle','docked'); hold on; grid on; box on;
scatter(Pin,TargetPout,sz,AGCErrorEval,'filled'); colorbar;
axis([Pin1-2 PoutMax-MinGain+2 Pin1+MaxGain-2 PoutMax+2])
xlabel('Pin (dBm)')
ylabel('Pout (dBm)')
title('AGC Error')
set(gca, 'FontSize',20)

figure('windowstyle','docked'); hold on; grid on; box on;
scatter(Pin,TargetPout,sz,NFEval,'filled'); colorbar;
axis([Pin1-2 PoutMax-MinGain+2 Pin1+MaxGain-2 PoutMax+2])
xlabel('Pin (dBm)')
ylabel('Pout (dBm)')
title('Noise Figure')
set(gca, 'FontSize',20)

