
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import numpy as np
import peakutils
import plotly.plotly as py
import plotly.graph_objs as go
import matplotlib.pyplot as plt
import copy
from scipy.signal import butter, lfilter, freqz


# this code works for sigle stage amps without scmd like: BOA4C211BDAH.

path_uut = '54_HOA2-HA-E21211-HA-2s/20190319_083641'  # it is the spectrum for the 2nd stage? 
mask_file = path_uut + '/mask/mask1.txt'
ampCh_0, srcCh_0, laser3Pwr_0, targetOutput_0, pin_0, targetGain_0 = [], [],[], [],[], []
 
# Import the files from mask.txt 
with open(mask_file) as f:
    for line in f:
        if(line.strip() != ''): 
            line_spl = line.split('\t') # vamos quebrar cada linha por tab
            cols = line_spl[0:]
                      
            targetGain_0.append(line_spl[0:1])
            pin_0.append(line_spl[1:2])
            targetOutput_0.append(line_spl[2:3])
            laser3Pwr_0.append(line_spl[3:4])
            srcCh_0.append(line_spl[4:5])
            ampCh_0.append(line_spl[5:6])
                        
# Function to change a list the data from string to float 
def numerizar (list):
  new_list = []
  #print(str(list))
  list.pop(0)  #delete the first element.  Delete the title of the first row.
  for i in range (0,len(list)):
      #print (list[i][0])
      new_list.append(float(list[i][0]))
  return new_list  


ampCh, srcCh, laser3Pwr, targetOutput, pin, targetGain = [], [],[], [],[], []

targetGain = numerizar(targetGain_0)    
ampCh = numerizar(ampCh_0)
pin = numerizar(pin_0)    

print (str(targetGain))
print (str(pin))


def ls(ruta = Path.cwd()):
    list = []
    #return [arch.name for arch in Path(ruta).iterdir() if arch.is_file()]
    for arch in Path(ruta).iterdir():
      if arch.is_file():
          name = arch.name
          list.append(name)   
    return list 


files = []
files = ls(path_uut + '/spectrum')


regis = path_uut.split('-')

# Number of stages... choose the current stage?
#print (regis[4][0:2])
stg = regis[4][0:2]

#print(files[0])
#print(files[2])

src_files = []
TargetGain = []
Pin = []
files_aux = []
TargetPout = []

for path in files:

    registro = path.split('_')
    if(registro[4]== "Amp.txt" ):
      TargetGain.append(float(registro[1]))
      Pin.append(float(registro[3]))
      TargetPout.append(float(registro[3])+float(registro[1]))
    else:
      src_files.append(path)   


for i in range(len(Pin)):
    files_aux.append([])

amp_files = []


# Re-organizing the columns as below:
#TargetGain| Pin | TargetPout
#
for i in range(len(Pin)):
    files_aux[i].append(TargetGain[i])
    files_aux[i].append(Pin[i])
    files_aux[i].append(TargetPout[i])

    amp_files.append(files_aux[i])

#print(amp_files)

# Sorting the values of power. 
valores_serie = pd.Series(amp_files)
Amp_files = valores_serie.sort_values(ascending=True)

#Re-defining the index of the list
Amp_files = Amp_files.reset_index(drop=True)

#print(Amp_files)

#print (len(Amp_files))

TargetGain_s = []
Pin_s = []
TargetPout_s = []

for file in Amp_files:
   TargetGain_s.append(file[0])
   Pin_s.append(file[1])
   TargetPout_s.append(file[2]) 

planck      = 6.626068e-34;
light       = 299792458;
res         = 0.097 * 1e-9; #0.1 * 1e-9

print ("res:  "+str(res))


# Define Fiber value according Stage.
if stg == "1s":
    #STG 1
    Fiber = 19.96
elif stg == "2s":
    #STG 2
    Fiber = 0




#print(max(lista))

#Exemplo ordenação
#-------------------
#mat = [[10,-17,-7],[10,-15,-6],[11,-18,-7],[10,-16,-6],[11,-15,-7]]
#valores_serie = pd.Series(mat)
#print(valores_serie)
#new=valores_serie.sort_values(ascending=True)
#print(new)
#-------------------


# maximum value of Gain and Pin of Src files since Amp spectrum files.  
file_spectrumSrc = path_uut+"/spectrum/Gain_%0.2f_Pin_%0.2f_Src.txt" % (max(TargetGain_s)-2, max(TargetPout) - max(TargetGain_s))
spectrumSrc_wave = []
spectrumSrc_power = []
print(file_spectrumSrc)

with open(file_spectrumSrc) as f:
    for line in f:
        if(line.strip() != ''): 
            line_spl = line.split('\t') # vamos quebrar cada linha por tab
            spectrumSrc_wave.append(line_spl[0:1])
            spectrumSrc_power.append(line_spl[1:2])


spectrumSrc_wave_n = []
spectrumSrc_power_n = []

#Change the data to float.
spectrumSrc_wave_n = numerizar(spectrumSrc_wave)
spectrumSrc_power_n = numerizar(spectrumSrc_power)

#Calculate the resolution
Dres = (np.max(spectrumSrc_wave_n) - np.min(spectrumSrc_wave_n)) / (len(spectrumSrc_wave_n)-1)
resol = np.ceil(3*res/Dres) #


#print(spectrumSrc_power_n)

print(spectrumSrc_power_n[845])

min = np.min(spectrumSrc_power_n)
max = np.max(spectrumSrc_power_n)
threshold = (min + max) / 2
#print("value max: " + str(max) + ", value min: " + str(min) +  ", threshold:  "+ str(threshold))
 

#print(len(peakutils.indexes(spectrumSrc_power_n, thres=threshold))) 

#for i in range(43011,43017): 
#  #print ("index: "+ str((i)*(1/100)))	
#  #print("tressss: " + str(len(peakutils.indexes(spectrumSrc_power_n, thres=(i+1)*(-1)))))
#  #picos = len(peakutils.indexes(spectrumSrc_power_n, thres=threshold, min_dist = -i, thres_abs = True ))
#  picos = len(peakutils.indexes(spectrumSrc_power_n, thres=(-i/1000), thres_abs = True, min_dist= max ))
# #print(tres)
#  #if(tres != 72 ):
#  if(picos != 0 ):
#     #print("diferente :" + str((i+1)*(1/100)) )
#     print("valueee: " + str(picos) + "  threshold= "+ str(-i/1000) + "listas: " + str(peakutils.indexes(spectrumSrc_power_n, thres=(-i/1000), thres_abs = True, min_dist= max))  + "\n") 


#print(len(peakutils.indexes(spectrumSrc_power_n, thres=)))


def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y


# Filter requirements.
order = 1
fs =  30      # sample rate, Hz
cutoff = 3  # desired cutoff frequency of the filter, Hz  # com 3 ok


# Filter the data, and plot both the original and filtered signals.
y = butter_lowpass_filter(spectrumSrc_power_n, cutoff, fs, order)


min = np.min(y)
max = np.max(y)
print( (min + max) / 2)

#listaa = peakutils.indexes(spectrumSrc_power_n, thres=-50, thres_abs = True, min_dist= max)

#list peaks positions of peaks.
listaa = peakutils.indexes(y, thres=((min + max) / 2), thres_abs = True, min_dist= max)

#print (str("longiiii "+ str(len(listaa)) + " lista filtrada indexes: " + str(listaa) ) + "\n")


AllLamS_0 = []
#AllLamS = spectrumSrc(loc,1)'*1e9
for i in range(len(listaa)):
    frec_max = spectrumSrc_wave_n[listaa[i]]
    
    #mult= float(frec_max)*1e9
    AllLamS_0.append(frec_max*1e9)
    

print(int(AllLamS_0[0]))
PoutPM = [];
PoutPD = [];
Src2PM_lst = [];
UUT2PM_lst = [];
flag = True
Spec_Src = []
Spec_Src_power = []

for i in range (0, len(TargetGain_s)):
    #pData = path_uut+"/result/Gain_%0.2f_Pin_%0.2f_PM.txt" % (max(TargetGain_s)-2, max(TargetPout) - max(TargetGain_s))
    if (Pin_s[i] < 0 or Pin_s[i] >9) :
      path_pData = path_uut+"/result/Gain_%0.2f_Pin_%0.2f_PM.txt" % (TargetGain_s[i], Pin_s[i])  #
    else:
      path_pData = path_uut+"/result/Gain_%0.2f_Pin_0%0.2f_PM.txt" % (TargetGain_s[i], Pin_s[i])  #

    with open(path_pData) as f:

          first_line = True
          for line in f:
             if(line.strip() != ''):
                 if(first_line == True ):
                    first_line = False
                 else: 
                    line_spl = line.split(' ') # vamos quebrar cada linha por espaco
                    #print( line_spl)
                    PoutPM.append(float(line_spl[1:2][0]))
                    
                    if(stg == "1s"):
                       PoutPD.append(float(line_spl[3:4][0])) 
                    elif(stg == "2s"):
                       PoutPD.append(float(line_spl[5:6][0])) 
                    
                    if(stg == "1s"):
                      #if( (line_spl[1:2]<-100) or (line_spl[8:9]<-45)):
                      if( (float(line_spl[1:2])>-101) or (float(line_spl[8:9])>-46)):
                          Src2PM_lst.append(float(line_spl[8:9][0])+Fiber-float(line_spl[0:1][0]))
                          UUT2PM_lst.append(float(line_spl[3:4][0])-float(line_spl[1:2][0]))
                    elif(stg == "2s"):
                      #if (line_spl[1:2][0]<-100):
                      if (float(line_spl[1:2][0])>-101):
                         Src2PM_lst.append(float(line_spl[4:5][0])-float(line_spl[0:1][0]))
                         UUT2PM_lst.append(float(line_spl[5:6][0])-float(line_spl[1:2][0]))
 
                    b = copy.deepcopy(TargetGain_s)
                    #maxx = max(TargetGain_s)
                    aux = pd.Series(b).max()
                    #print("maxxx "+ str(aux))

                    if( (TargetGain_s[i]== pd.Series(TargetGain_s).max()-8) and (TargetPout_s[i] > pd.Series(TargetPout_s).max()-9) and flag == True ):
                      flag = False
                      file_spectrumSrc = path_uut+"/spectrum/Gain_%0.2f_Pin_%0.2f_Src.txt" % (TargetGain_s[i], Pin_s[i])
                      with open(file_spectrumSrc) as f:
                          first_line = True
                          for line in f:
                              if(first_line == True ):
                                 first_line = False
                              else:
                                 if(line.strip() != ''): 
                                     line_spl = line.split('\t') # vamos quebrar cada linha por tab
                                     Spec_Src_power.append(line_spl[1:2][0])

factor = 0.000215*len(Spec_Src_power)/0.1


#conversion of dBm to watt.  It receive a list of numbers and return a list of numbers in Watt.
def conv_dbm2w(list_db):
  expo_var = [] 
  src_power = []
  for element in list_db:
       src_power.append(float(element)/10) 
  expo_var = (np.power(10, (src_power)) * 1e-3) 
  return expo_var


# reproducing eq. power  = 10*log10(trapz(10.^(spectrumSrc(:,2)/10)*1e-3) *1/factor*1e3); 
#calculating integral
expo_var2 = np.trapz(conv_dbm2w(Spec_Src_power) )  #ok
#print("integral  "+ str(expo_var2))

power = 10*np.log10(expo_var2 *(1/factor)*1e3) 
print("integral  "+ str(power))


OSA2PM = float(line_spl[0:1][0]) - power
print("OSA2PM: "+ str(OSA2PM))


#print(PoutPD)
#print(Src2PM)
print(len(Spec_Src))
media = pd.Series(Src2PM_lst)
Src2PM = media.mean()
print ("Src2PM: "+ str(Src2PM))

media = pd.Series(UUT2PM_lst)
UUT2PM = media.mean()
print ("UUT2PM: "+ str(UUT2PM))


print (len(TargetGain_s))
print (len(Spec_Src_power))
Src2UUT = 0;
list_list_src_power = []


# this function extract datas from files in spectrum folder.
def import_spec(data_path):
   with open(data_path) as f:
       spec_list_power = []
       spec_list_frec = []
       first_line = True
       for line in f:
           if(first_line == True ):
              first_line = False
           else:
              if(line.strip() != ''): 
                  line_spl = line.split('\t') # vamos quebrar cada linha por tab
                  spec_list_power.append(float(line_spl[1:2][0])) 
                  spec_list_frec.append(float(line_spl[0:1][0])) 
   #it returns 2 columns frecuency and power
   return spec_list_frec, spec_list_power


numbers = []  
spec_src_power_all = []
spec_amp_power_all = [] 

spec_amp_frec_all = [] 
spec_src_frec_all = [] 

for i in range (0, len(TargetGain_s)):  # 550 files 
   #spectrumAmp = dlmread([[path '\spectrum\'] sprintf('Gain_%05.2f_Pin_%05.2f_Amp.txt',TargetGain(j),Pin(j))],'	',2,0);
   if (Pin_s[i] < 0 or Pin_s[i] >9) :
       #file_spectrumSrc = path_uut+"/spectrum/Gain_%0.2f_Pin_%0.2f_Src.txt" % (TargetGain_s[i], TargetPout[i] -TargetGain_s[i])
       file_spectrumSrc = path_uut+"/spectrum/Gain_%0.2f_Pin_%0.2f_Src.txt" % (TargetGain_s[i], Pin_s[i])
       file_spectrumAmp = path_uut+"/spectrum/Gain_%0.2f_Pin_%0.2f_Amp.txt" % (TargetGain_s[i], Pin_s[i])
   else:
       #file_spectrumSrc = path_uut+"/spectrum/Gain_%0.2f_Pin_0%0.2f_Src.txt" % (TargetGain_s[i], TargetPout[i] -TargetGain_s[i])
       file_spectrumSrc = path_uut+"/spectrum/Gain_%0.2f_Pin_0%0.2f_Src.txt" % (TargetGain_s[i], Pin_s[i])
       file_spectrumAmp = path_uut+"/spectrum/Gain_%0.2f_Pin_0%0.2f_Amp.txt" % (TargetGain_s[i], Pin_s[i])
   numbers.append(str(TargetGain_s[i]) + " " + str(Pin_s[i]))

   spectrum_src_power = []   
   spectrum_amp_power = []   
 
   spec_src_aux = import_spec(file_spectrumSrc)  
   for element in spec_src_aux[1]: 
      spectrum_src_power.append(element +  (Src2PM - Src2UUT + OSA2PM)) 

   #print ("factor : "+ str(factor0 ))
   spec_amp_aux = import_spec(file_spectrumAmp)  
   for element in spec_amp_aux[1]: 
      spectrum_amp_power.append(element +  (UUT2PM + OSA2PM))

   #listas de debug
   spec_src_power_all.append(spectrum_src_power) # just power.  All the 550 files. 
   spec_amp_power_all.append(spectrum_amp_power) 
   
   spec_src_frec_all.append(spec_src_aux[0]) 
   spec_amp_frec_all.append(spec_amp_aux[0]) 

   #sltRange = spectrumAmp(:,1)*1e9>min(pmaxSrc(:,1))-7 & spectrumAmp(:,1)*1e9<max(pmaxSrc(:,1))+7;  NEXT vector of "1s" 


# Filter the data, and plot both the original and filtered signals.
# It calculates over the first src file. Just for test. 
y = butter_lowpass_filter(spec_src_power_all[0], cutoff, fs, order)

pmaxSrc_power = []
pmaxSrc_frec = []

min = np.min(y)
max = np.max(y)
#print( (min + max)*0.7 )

loc = peakutils.indexes(y, thres=((min + max)*0.22 ), thres_abs = True, min_dist= max)
print ("looooc" + str(loc))


#it adds to the list the power and frec of peaks
for elem in loc:
  pmaxSrc_power.append( spec_src_power_all[0][elem] )
  pmaxSrc_frec.append( spec_src_frec_all[0][elem]*1e9 )

#it calculates the power and frec of peakes   
print (pmaxSrc_power)    
print (pmaxSrc_frec)    

#sltRange = spectrumAmp(:,1)*1e9>min(pmaxSrc(:,1))-7 & spectrumAmp(:,1)*1e9<max(pmaxSrc(:,1))+7;
new_list = [] 
filtered_list = [] 
sltSpectrumAmp_power = []
sltSpectrumAmp_frec = [] 

for i in range (0, len(spec_amp_power_all[0])):  # in the file  
   
   if (spec_amp_frec_all[0][i]*1e9 > (np.min(pmaxSrc_frec)-7) and spec_amp_frec_all[0][i]*1e9 < (np.max(pmaxSrc_frec)+7 ))  :
       new_list.append(1)
       filtered_list.append(i)
       sltSpectrumAmp_power.append(spec_amp_power_all[0][i] )
       sltSpectrumAmp_frec.append(spec_amp_frec_all[0][i] ) 
 
   else:  
       new_list.append(0) 
       #print (spec_amp_frec_all[0][i]*1e9) 
       #print (np.min(pmaxSrc_frec)) 
#print ("newww" + str(new_list) )
#print ("newww" + str(filtered_list) )


frec = []
deltafrec = []

 
#DeltaFreq   = (Freq.^2).*res./light;
for elem in pmaxSrc_frec:
   frec.append( (light * 1e9)/elem)
   deltafrec.append( ( np.power( ((light * 1e9)/elem),2) * res / light ) ) 

print("frec: "+ str(frec) )
print("deltafrec: "+ str(deltafrec) )

factor1  = 0.000215 * len(spec_amp_power_all[0])/(0.1)

print ("factor : "+ str(factor1 ))

sltSpectrumSrc = [] 
sltSpectrumAmp = [] 

#sltSpectrumSrc = spectrumSrc(sltRange,1);
for elem in filtered_list:
     sltSpectrumSrc.append(spec_src_power_all[0][elem]) # este aqui para q????? 

print ("sltSpectrumSrc***" + str( len( sltSpectrumSrc  )) )
#print ("*************" + str( spec_src_power_aux_1[1]) )
#print("lista ===  "+ str(numbers))
 

#  reproducing : id = find(sltSpectrumAmp(:,2)==min(sltSpectrumAmp(locMax(i):locMax(i+1),2)));
# Filter the data, and plot both the original and filtered signals.
# It calculates over the first src file. Just for test. 
#yy = butter_lowpass_filter(sltSpectrumSrc, cutoff, fs, order)
yy = butter_lowpass_filter(sltSpectrumAmp_power, cutoff, fs, order)
 
min = np.min(y)
max = np.max(y)

locMax = peakutils.indexes(yy, thres=((min + max)*0.22 ), thres_abs = True, min_dist= max)
print( "locMax: " + str(locMax))

minim = np.min(sltSpectrumAmp_power[locMax[0]:locMax[1]]) 
print( "minim: " + str(minim))

#1. next step procurr pontos minimos pocisões  entre esse range total de sltspectrumAMP. 
#so para dois maximos, fazer para 3 maximos ou mais!
 
loc_0 = []
for i in range(locMax[0],locMax[1]):
  if(sltSpectrumAmp_power[i] == minim ):
     loc_0.append(i)

pmaxAmp_pwr = [] 
pmaxAmp_frec = [] 
for elem in locMax:
   pmaxAmp_pwr.append(sltSpectrumAmp_power[elem])
   pmaxAmp_frec.append(sltSpectrumAmp_frec[elem]*1e9)
   #pmaxAmp_pwr.append(yy[elem])
print ("pmaxAMP :" + str(pmaxAmp_pwr) )

#Calculate the resolution again. Maybe optional?
#Dres_0 = (np.max(spec_src_frec_all[0]) - np.min(spec_src_frec_all[0])) / (len(spec_src_frec_all[0])-1)
#resol = np.ceil(3*res/Dres_0) # 7
#print("resl = "+ str(resol)  )

#  ---  extremo esquerdo  ---- 
#id = find(sltSpectrumAmp(:,2)==min(sltSpectrumAmp(locMax(1)-ceil(5*res/Dres):locMax(1),2)));
minim_izq  = np.min(sltSpectrumAmp_power[(locMax[0]-int(resol)):locMax[0]])
#print ("minn" + str(minim_izq) )
#id(id<(locMax(1)-ceil(5*res/Dres)) | id>locMax(1)) = [];
for i in range(locMax[0]-int(resol),locMax[0]):
  if(sltSpectrumAmp_power[i] == minim_izq ):

     loc_0.append(i)


#  ---  extremo direto  ---- 
#id = find(sltSpectrumAmp(:,2)==min(sltSpectrumAmp(locMax(end):locMax(end)+ceil(5*res/Dres),2)));
minim_der =  np.min(sltSpectrumAmp_power[(locMax[len(locMax)-1]):( locMax[len(locMax)- 1]  + int(resol) ) ])
#id(id<locMax(end) | id>(locMax(end)+ceil(5*res/Dres))) = [];
for i in range((locMax[len(locMax) - 1]) , (locMax[len(locMax) -1 ]+ int(resol) ) ):
  if(sltSpectrumAmp_power[i] == minim_der ):
     loc_0.append(i)
print("loc_0: " + str(loc_0)  )

# pminAmp = [sltSpectrumAmp(loc,1)*1e9 sltSpectrumAmp(loc,2)];
pminAmp_pwr = []
for elemen in loc_0:
   pminAmp_pwr.append ( sltSpectrumAmp_power[elemen])   

print("pminAmp= " + str(pminAmp_pwr) )

PaseAuxIN = []
pminSrc_pwr = []
#
#pminSrc = [sltSpectrumSrc(loc,1)*1e9 sltSpectrumSrc(loc,2)];


for elem in loc_0:
   pminSrc_pwr.append(sltSpectrumSrc[elem] )
print("pminSrc_pwr" + str(pminSrc_pwr ) )
  
#  
#PaseAuxIN = ([10.^((pminSrc(:,2)+NoiseConv)/10);0] + [0;10.^((pminSrc(:,2)+NoiseConv)/10)])./2;  # esta muy raro esto.  perguntar o Joao

first_p = conv_dbm2w(pminSrc_pwr)
var1 = []
for elem in first_p:
  var1.append(float(elem)*1e3) 
var1.append(0) 
print ("var1: " + str(var1)  )

var2 = np.roll(var1,1)

print ("var2 ;" + str(var2)  )
PaseAuxIN = []
for elem in range (0, len(var2)):
   # calculates the mean
   PaseAuxIN.append ( (var1[elem] + var2[elem])/2  ) 

print (PaseAuxIN) 
	
# lamS = pmaxAmp(:,1);
lamS = []
for elem in pmaxAmp_frec:
   lamS.append(elem)
#print( str(lamS) ) 


#***        PaseIN  = [lamS PaseAuxIN(2:end-1)*1e-3];  
aux_1 = []

for i in range(1, len(PaseAuxIN)-1 ): #Paseaux sempre vai ter 4 elem?
   aux_1.append(PaseAuxIN[i]*1e-3 ) 

PaseIN = []
PaseIN.append(lamS)
PaseIN.append(aux_1) 
print (PaseIN)




# --****        PaseAuxOUT= ([10.^((pminAmp(:,2)+NoiseConv)/10);0] + [0;10.^((pminAmp(:,2)+NoiseConv)/10)])./2;
first_p = []
first_p = conv_dbm2w(pminAmp_pwr)
var1 = []
for elem in first_p:
  var1.append(float(elem)*1e3) 
var1.append(0) 
#print ("var1: " + str(var1)  )
var2 = np.roll(var1,1)
#print ("var2 ;" + str(var2)  )

PaseAuxOUT = []
for elem in range (0, len(var2)):
   # calculates the mean
   PaseAuxOUT.append ( (var1[elem] + var2[elem])/2  ) 


# --****       PaseOUT = [lamS PaseAuxOUT(2:end-1)*1e-3];
PaseOUT = []
aux_1 = []

for i in range(1, len(PaseAuxOUT)-1 ): #Paseaux sempre vai ter 4 elem?
   aux_1.append(PaseAuxOUT[i]*1e-3 ) 


PaseOUT.append(lamS)
PaseOUT.append(aux_1) 
print (PaseOUT)




#for i in range(4100,5300): 
#  #print ("index: "+ str((i)*(1/100)))	
#  #print("tressss: " + str(len(peakutils.indexes(spectrumSrc_power_n, thres=(i+1)*(-1)))))
#  #picos = len(peakutils.indexes(spectrumSrc_power_n, thres=threshold, min_dist = -i, thres_abs = True ))
#  picos = len(peakutils.indexes(y, thres=(-i/100), thres_abs = True, min_dist= max ))
#  #print(tres)
#  #if(tres != 72 ):
#  if(picos != 0 ):
#     #print("diferente :" + str((i+1)*(1/100)) )
#     print("valueee: " + str(picos) + "  threshold= "+ str(-i/100))



#for i in range(100):
	#print(i)
#    print(peakutils.indexes(spectrumSrc_power_n, thres=(np.min(spectrumSrc_power_n)-np.max(spectrumSrc_power_n))/2)

#trace = go.Scatter(
#    x = [j for j in range(len(spectrumSrc_power_n))],
#    y = spectrumSrc_power_n,
#    mode = 'lines'
#)

#data = [trace]
#py.iplot(data, filename='milk-production-plot')

# gca stands for 'get current axis'
#ax = plt.gca()


#for i, row in df.iterrows():

#plt.plot(row['list1'], row['list2'])

#df.plot(kind='line',x='name',y='num_children',ax=ax)




ind = np.arange(len(spectrumSrc_power_n))

fig = plt.figure()
#fig_f = plt.figure()

ax = fig.add_subplot(111)
ax.plot(ind, spectrumSrc_power_n)

ax1 = fig.add_subplot(111)
ax1.plot(ind, y)

#for i in range(len(listaa)):
#   ax1.plot( listaa[i], spectrumSrc_power_n[listaa[i]], 'go') # green bolinha

varrr = 'exemm'	
plt.savefig(varrr+".svg")

#print(spectrumSrc_wave_n)
#print(spectrumSrc_power_n)

