#Start with import statements 
import pandas as pd
import glob

#Read in all variable forecast files
path = '/Users/ellendyer/Desktop/Nov23_Python/collaborative_python/'+'new_files/'
stpr = pd.read_excel(glob.glob(path+"Awash*R*.xlsx")[0],header=None,index_col=False)
stta = pd.read_excel(glob.glob(path+"Awash*Max*.xlsx")[0],header=None,index_col=False)
stti = pd.read_excel(glob.glob(path+"Awash*Min*.xlsx")[0],header=None,index_col=False)

#Check to see if there is the WRF header row with one column
#and fix column headers
if stpr[0][0].startswith('WRF') == True:
    stpr = stpr.drop(0)
cols = stpr.iloc[0]
stpr = stpr[1:]
stpr.columns = cols
if stta[0][0].startswith('WRF') == True:
    stta = stta.drop(0)
cols = stta.iloc[0]
stta = stta[1:]
stta.columns = cols
if stti[0][0].startswith('WRF') == True:
    stti = stti.drop(0)
cols = stti.iloc[0]
stti = stti[1:]
stti.columns = cols


coln = len(stpr.columns)
for S in stpr.STATION:
    sinfo_pr = stpr[stpr['STATION'] == S]
    sinfo_ta = stta[stta['STATION'] == S]
    sinfo_ti = stti[stti['STATION'] == S]
    
    print('** At ',S,' **')
    print('The rainfall forecast is:')
    for C in range(3,coln-1):
        print('. ',sinfo_pr.columns[C], "{:.2f}".format(sinfo_pr[sinfo_pr.columns[C]].values[0]), 'mm')
    print('. ',sinfo_pr.columns[-1], 'over this period:', "{:.2f}".format(sinfo_pr[sinfo_pr.columns[-1]].values[0]), 'mm')
    
    print('The maximum temperature forecast is:')
    for C in range(3,coln-1):
        print('. ',sinfo_ta.columns[C], sinfo_ta[sinfo_ta.columns[C]].values[0], 'C')
    
    print('The minimum temperature forecast is:')
    for C in range(3,coln-1):
        print('. ',sinfo_ti.columns[C], sinfo_ti[sinfo_ti.columns[C]].values[0], 'C')
        
    print(" ")
    print("---------------------------------------------")
    print(" ")
    
    

