from pathlib import Path
import pandas as pd
import numpy
#from scipy.constants import pi


path = '54_HOA2-HA-E21211-HA-2s/20190319_083641'
mask_file = path + '/mask/mask1.txt'

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
files = ls(path + '/spectrum')


regis = path.split('-')

# Number of stages
#print (regis[4][0:2])
stg = regis[4][0:2]



print(files[0])
print(files[2])

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



#IL
if stg == "1s":
    #STG 1
    Fiber = 19.96
elif stg == "2s":
    #STG 2
    Fiber = 0


#Exemplo ordenação
#-------------------
#mat = [[10,-17,-7],[10,-15,-6],[11,-18,-7],[10,-16,-6],[11,-15,-7]]
#valores_serie = pd.Series(mat)
#print(valores_serie)
#new=valores_serie.sort_values(ascending=True)
#print(new)
#-------------------

NoiseConv = 0






