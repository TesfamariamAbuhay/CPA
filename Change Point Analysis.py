# -*- coding: utf-8 -*-
"""
Created on Wed May 23 12:14:40 2018

@author: user
"""
import matplotlib.pylab as plt
import matplotlib
import pandas as pd
import numpy as np
import seaborn as sns
Tropon = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Troponin.csv")
Data = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Clinical Test.csv")
ACS_Data = Data[['Department_event', 'sex', 'patient_id','date_time_test', 'result', 'condition','tp', 'tropon','crea', 'chlor']]
ACS_Data['Date'] = pd.to_datetime(ACS_Data.date_time_test)
ACS_Data.index = ACS_Data.Date
ACS_Data_ID = ACS_Data.patient_id.value_counts()
columns = ['Date', 'Troponin'] #'Creatinin'
ID = ACS_Data_ID[ACS_Data_ID > 30]
Troponin_Creatnin1 = pd.DataFrame(columns = columns)
for i in ID.index:
    Case2 = ACS_Data[ACS_Data.patient_id == i]
    Case2['Date'] = pd.to_datetime(Case2.date_time_test)
    Case2.index = Case2.Date
    Case2_Creat = Case2[pd.to_numeric(Case2['crea'], errors='coerce').notnull()]
    Case2_Trop = Case2[pd.to_numeric(Case2['tropon'], errors='coerce').notnull()]
    Case2_Trop.tropon = Case2_Trop.tropon.astype('float')
    Case2_Trop = Case2_Trop[Case2_Trop.tropon > 0.01]
    Case2_Trop = Case2_Trop[Case2_Trop.tropon < 50]
    Case2_Trop = Case2_Trop.tropon
    #Case2_Trop = Case2_Trop.tropon.resample('B', how = 'mean').dropna()
    Case2_Creat.crea = Case2_Creat.crea.astype('float')
    Case2_Creat = Case2_Creat[Case2_Creat.crea > 18.0]
    Case2_Creat = Case2_Creat[Case2_Creat.crea < 800]
    Case2_Creat = Case2_Creat.crea
    #Case2_Creat = Case2_Creat.crea.resample('B', how = 'mean').dropna()
    #cut = np.min([len(Case2_Trop), len(Case2_Creat)])
    #if cut == 0:
    #    pass
    #else:
    #Case2_Creat2 = Case2_Creat[Case2_Creat.index.isin(Case2_Trop.index)]
    #Case2_Trop2 = Case2_Trop[Case2_Trop.index.isin(Case2_Creat.index)]
    #Troponin_Creatnin = pd.DataFrame([list(Case2_Creat2.index), list(Case2_Trop2.values)[:len(Case2_Trop2)], list(Case2_Creat2.values)[:len(Case2_Creat2)]]).T
    if len(Case2_Trop) > 5:
    #Case2_Trop2 = Case2_Trop[Case2_Trop.index.isin(Case2_Creat.index)]
        Troponin_Creatnin = pd.DataFrame([list(Case2_Trop.index), list(Case2_Trop.values)[:len(Case2_Trop)]]).T
        Troponin_Creatnin.columns = columns
        Troponin_Creatnin['Patient_ID'] = i
        Troponin_Creatnin1= Troponin_Creatnin1.append(Troponin_Creatnin)
    else:
        pass
    

Troponin_Creatnin1.to_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Creatinin_Only.csv")
Tropon = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Troponin_Only.csv")
Patient_Info = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Patient_Info_Age_Sex.csv")
Patient_Info.Death = Patient_Info.Death.replace([0,1], ['Alive', 'Died'])
#Tropon = Tropon.drop('Date.1', axis = 1)
#ID = Tropon.Patient_ID.value_counts()
#ID = ID[ID > 7]
#CUSUM for Troponin only
columns_C = ['Date','Patient_ID',  'Troponin', 'CUSUM', 'Age', 'Gender', 'Death']
CU_SUM_Mid_T = pd.DataFrame(columns = columns_C)
Mean_Mid = []
CU_SUM = []
t = 0
for c in Tropon.Patient_ID.unique():
    Data = Tropon[Tropon.Patient_ID == c].dropna()
    try:
        patient_info = list(Patient_Info[Patient_Info.patient_id == c][['age', 'sex', 'Death']].values[0])
    except IndexError:
        patient_info = [67, 'Male', 0]
    Data.index = Data.Date
    mean = Data.Troponin.dropna().mean()
    Mean_Mid.append(mean)
    CU_SUM.append([0])
    res = pd.DataFrame([Data.index[0], c, 0, 0, patient_info[0], patient_info[1], patient_info[2]]).T
    res.columns = columns_C
    CU_SUM_Mid_T = CU_SUM_Mid_T.append(res)
    d = 0
    for i in Data['Troponin'].values:
        s = CU_SUM[t][-1] + (i - mean)
        CU_SUM[t].append(s)
        res = pd.DataFrame([Data.index[d], c, i, s, patient_info[0], patient_info[1], patient_info[2]]).T
        res.columns = columns_C
        CU_SUM_Mid_T = CU_SUM_Mid_T.append(res)
        d = d + 1
    t = t + 1
#CU_SUM_Mid_T.index = CU_SUM_Mid_T.Date
CU_SUM_Mid_T['Date']= pd.to_datetime(CU_SUM_Mid_T.Date)
CUSUM = pd.DataFrame()
for i in CU_SUM_Mid_T.Patient_ID.unique():
    Data = CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == i]
    shift = len(Data)
    Data['Diff_Min'] = ((Data['Date'] - Data['Date'].shift(1)).dt.seconds)/60
    Data['Inter-Test Hour'] = Data['Diff_Min']/60
    LOS  = []
    Data = Data.fillna(0)
    for t in range(shift):
        LOS.append(np.sum(Data['Inter-Test Hour'][:t + 1]))
    Data['LoS in Hours'] = LOS
    CUSUM = CUSUM.append(Data)
CUSUM = CUSUM.fillna(0)
CUSUM.index = CUSUM['LoS in Hours']
CUSUM.Death = CUSUM.Death.replace([0,1], ['Discharged', 'Died'])
CU_SUM_Mid_T.to_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Troponin_CUSUM.csv")

plt.figure(figsize = (14, 8)), 
CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == 432895].CUSUM.plot(), 
label = [round(i) for i in list(CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == 432895].LoS.values)]
xticks = CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == 432895].LoS.index
plt.xticks(xticks, label)

#age category
age = []
for i in CU_SUM_Mid_T.Age:
    if i < 55:
        age.append('<55')
    elif i in range(55,65):
        age.append('55-64')
    elif i in range(65, 75):
        age.append('65-75')
    else:
        age.append('>75')
CU_SUM_Mid_T['Age_Group'] = age

#plotting troponin distribution in age and gender categories
import matplotlib as mpl
mpl.style.use('classic')
plt.figure(figsize=(14, 8))
sns.boxplot(x = 'Age_Group', y = 'Troponin', hue = 'Gender', data = CU_SUM_Mid_T)
plt.figure(figsize=(14, 8))
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Age Group', fontsize=20)
plt.ylabel('Troponin', fontsize=20)
plt.legend(fontsize = 20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Age Group', fontsize=20)
plt.ylabel('Inter-Test Hour', fontsize=20)
sns.factorplot(x="Gender", y='Inter-Test Hour', col="Age_Group", data=CU_SUM_Mid_T, kind="violin", size=4, aspect=.5, col_order = col_order);
plt.ylim(0,)
plt.savefig('D:\\Simulation Model\\Change Point Analysis\\Data\\factorplot_inter-test hour.png', format='png', dpi=300, bbox_inches="tight")


ax = plt.subplot(121), 
plt.plot(x = CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == 383832]['Inter-Test Hour'], y = CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == 383832].CUSUM), 
t = ax.get_yticks()

#Troponin at Admission
Tro_Adm = [CU_SUM_Mid_T[CU_SUM_Mid_T.Patient_ID == i].head(3).values[1] for i in CU_SUM_Mid_T.Patient_ID.unique()]
Tro_Adm = pd.DataFrame(Tro_Adm, columns = CU_SUM_Mid_T.columns)
col_order = ['<55', '55-64', '65-75', '>75']
plt.figure(figsize=(14, 8))
sns.boxplot(x = 'Age_Group', y = 'Troponin', hue = 'Gender', data = Tro_Adm, col_order = col_order)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Age Group', fontsize=20)
plt.ylabel('Troponin', fontsize=20)
plt.legend(fontsize = 20)
plt.xticks(fontsize=20)
plt.yticks( fontsize=20)
g = sns.factorplot(x="Gender", y='Troponin', col="Age_Group", data=Tro_Adm, kind="box", size=6, aspect=.6, col_order = col_order);
g.set_axis_labels("Gender", "Troponin");
g.fig.subplots_adjust(wspace=.02, hspace=.02);
plt.savefig('D:\\Simulation Model\\Change Point Analysis\\Data\\factorplot.png', format='png', dpi=300, bbox_inches="tight")

#hist of Troponin in age and gender category
plt.figure(figsize=(14, 8))
g = sns.FacetGrid(Tro_Adm, col="Age_Group", hue="Gender", col_wrap=2, col_order = col_order)
g.map(plt.hist, "Troponin")
g.add_legend(fontsize=20);
g.xticks(fontsize=20)
g.yticks( fontsize=20)

#Sample CUSUM Chart for both dischared and died patients
Patients = Tropon.Patient_ID.value_counts()
plt.figure(figsize = (14, 8))
fig = ['A)', 'B)', 'C)', 'D)', 'E)', 'F)', 'G)']
c = 231
f = 0
for i in Patients[:6].index:
    ax = plt.subplot(c) 
    CUSUM[CUSUM.Patient_ID == i].CUSUM.plot(marker = '<', ms = 10, color = 'green', lw = 3), 
    pp = list(Patient_Info[Patient_Info.patient_id == i][['age', 'sex', 'Death']].values[0])    
    plt.title('Age=' + str(pp[0]) + ', Gender=' + pp[1], fontsize = 20)
    plt.xticks(fontsize = 20, rotation = 'vertical'), plt.yticks(fontsize = 20), 
    if c == 231 or c == 232 or c == 233:
        ax.set_xlabel('')
    else:
        plt.xlabel('LoS in Hours', fontsize = 20)
    if c == 231 or c == 234:
        plt.ylabel('Troponin', fontsize = 20)
    else:
        ax.set_ylabel('')
    c = c + 1
    x = ax.get_xticks()
    y = ax.get_yticks()
    #plt.text(x[-3] + 1, y[-2], fig[f], color = 'red', fontsize = 20)
    f = f + 1
    #plt.locator_params(numticks=4)
plt.tight_layout()
plt.savefig('D:\\Simulation Model\\Change Point Analysis\\Data\\CUSUM Sample2.png', format='png', dpi=300, bbox_inches="tight")



#plotting CUSUM died patients
Died = CUSUM[CUSUM.Death == 1].Patient_ID.unique()
plt.figure(figsize = (14, 8))
ax = plt.subplot(231) 
CUSUM[CUSUM.Patient_ID == Died[0]].CUSUM.plot(marker = '*', ms = 10), 
plt.title(str(list(Patient_Info[Patient_Info.patient_id == Died[0]][['age', 'sex', 'Death']].values[0][:2])), fontsize = 20)
plt.xticks(fontsize = 20), plt.yticks(fontsize = 20), ax.set_xlabel('')
ax = plt.subplot(232)
CUSUM[CUSUM.Patient_ID == Died[4]].CUSUM.plot(marker = '*', ms = 10),
plt.title(str(list(Patient_Info[Patient_Info.patient_id == Died[4]][['age', 'sex', 'Death']].values[0][:2])), fontsize = 20) 
plt.xticks(fontsize = 20), plt.yticks(fontsize = 20), ax.set_xlabel('')
ax = plt.subplot(233)
CUSUM[CUSUM.Patient_ID == Died[-3]].CUSUM.plot(marker = '*', ms = 10)
plt.title(str(list(Patient_Info[Patient_Info.patient_id == Died[-3]][['age', 'sex', 'Death']].values[0][:2])), fontsize = 20)
plt.xticks(fontsize = 20), plt.yticks(fontsize = 20), ax.set_xlabel('')
plt.subplot(234),CUSUM[CUSUM.Patient_ID == Died[1]].CUSUM.plot(marker = '*', ms = 10)
plt.title(str(list(Patient_Info[Patient_Info.patient_id == Died[1]][['age', 'sex', 'Death']].values[0][:2])), fontsize = 20)
plt.xticks(fontsize = 20), plt.yticks(fontsize = 20), plt.xlabel('LoS in Hours', fontsize = 20)
plt.subplot(235),CUSUM[CUSUM.Patient_ID == Died[10]].CUSUM.plot(marker = '*', ms = 10)
plt.title(str(list(Patient_Info[Patient_Info.patient_id == Died[10]][['age', 'sex', 'Death']].values[0][:2])), fontsize = 20)
plt.xticks(fontsize = 20), plt.yticks(fontsize = 20), plt.xlabel('LoS in Hours', fontsize = 20)
plt.subplot(236),CUSUM[CUSUM.Patient_ID == Died[14]].CUSUM.plot(marker = '*', ms = 10)
plt.title(str(list(Patient_Info[Patient_Info.patient_id == Died[14]][['age', 'sex', 'Death']].values[0][:2])), fontsize = 20)
plt.xticks(fontsize = 20), plt.yticks(fontsize = 20),plt.xlabel('LoS in Hours', fontsize = 20)
plt.tight_layout()
plt.savefig('D:\\Simulation Model\\Change Point Analysis\\Data\\CUSUM for Died Patients.png', format='png', dpi=300, bbox_inches="tight")



#Calculating confidenc interval and level
#calculte estimation magnitude
estimator_magnitude_change_mid = []
for p in CUSUM.Patient_ID.unique():
    Case = CUSUM[CUSUM.Patient_ID == p]
    estimator_magnitude_change_mid.append(Case.CUSUM.max() - Case.CUSUM.min())
len(estimator_magnitude_change_mid)

#bootstraping
bootstrap_mid = []
t = 0
for i in Tropon.Patient_ID.unique():
    Case = Tropon[Tropon.Patient_ID == i]
    bootstrap_mid.append([])
    for b in range(1000):
        bootstrap_mid[t].append(np.random.choice(Case.Troponin, size = len(Case.Troponin)))
    t = t + 1
len(bootstrap_mid[0][0])

Boot_CU_SUM_Mid = []
Boot_Mean_Mid = []
t = 0
for c in range(len(bootstrap_mid)):
    Boot_CU_SUM_Mid.append([])
    Boot_Mean_Mid.append([])
    for b in range(1000):
        mean = np.mean(bootstrap_mid[c][b])
        Boot_Mean_Mid[t].append(mean)
        Boot_CU_SUM_Mid[t].append([0])
        for i in range(len(bootstrap_mid[c][b])):
            s = Boot_CU_SUM_Mid[t][b][-1] + (bootstrap_mid[c][b][i] - mean)
            Boot_CU_SUM_Mid[t][b].append(s)
    t = t + 1

boot_estimator_magnitude_change_mid = []
conf_inter_mid = []
diff_1000 = []
t = 0
for i in range(len(Boot_CU_SUM_Mid)):
    conf_inter_mid.append([])
    boot_estimator_magnitude_change_mid.append([])
    diff_1000.append([])
    for b in range(1000):
        #conf_inter_mid[i].append([])
        #boot_estimator_magnitude_change_mid[i].append([])
        #for topic in Boot_CU_SUM_Mid[i][b]:
        topic = Boot_CU_SUM_Mid[i][b]
        boot_dif_mid = np.max(topic) - np.min(topic)
        boot_estimator_magnitude_change_mid[i].append(boot_dif_mid)
        diff_1000[i].append(boot_dif_mid)
        conf_inter_mid[i].append(boot_dif_mid < estimator_magnitude_change_mid[i])
    t = t + 1

# plottting bootsraps and original data
plt.figure(figsize = (14, 8)) 
ax = plt.subplot()
plt.plot(CUSUM[CUSUM.Patient_ID == 432895].CUSUM.values, label = 'Original order', marker = markers[-1], lw = 3, ms =10),
[plt.plot(Boot_CU_SUM_Mid[0][i], label = 'Bootstrap' + str(i+1), marker = markers[i], lw = 3, ms =10) for i in range(6)], 
plt.legend(fontsize = 20), plt.xlim(0,55), 
plt.ylabel('CUSUM', fontsize = 20),
plt.xlabel('LoS in Hours', fontsize = 20),
x = ax.get_xticks() 
xt = [0, 100, 200, 250, 300, 400, 450]
plt.xticks(x, xt, fontsize = 20), 
plt.yticks(fontsize = 20)
plt.savefig('D:\\Simulation Model\\Change Point Analysis\\Data\\CUSUM bootstrap.png', format='png', dpi=300, bbox_inches="tight")


#Conf_Inter = pd.DataFrame(conf_inter_mid)
Id = []
change = []
t = 0
for i in CUSUM.Patient_ID.unique():
    Case = CUSUM[CUSUM.Patient_ID == i]
    CI = pd.DataFrame(conf_inter_mid[t], columns = ['Bootstrap' + str(i)])
    res = CI[CI.columns[0]].value_counts()
    if res.index[0] == 'False':
        interval = CI[CI.columns[0]].value_counts()[0]/1000*100
        #change.append(interval)
    else:
        interval = CI[CI.columns[0]].value_counts()[1]/1000*100
        #change.append(interval)
        
    if interval > 90:
        m = list(abs(Case.CUSUM))
        s = np.sort(m)
        change.append([i, np.max(m), m.index(np.max(m)),  m.index(s[-2]), m.index(s[-3]), round(interval)])
        #pp = list(Patient_Info[Patient_Info.patient_id == i][['age', 'sex', 'Death']].values[0])    
        #print('Age=' + str(pp[0]) + ', Gender=' + pp[1] + ', Result=' + pp[2])        
        Id.append(i)
    else:
        pass
    t = t + 1

#Change
plt.figure(figsize = (14, 8))
c = 231
f = 0
for i in Id[:6]:
    ax = plt.subplot(c) 
    CUSUM[CUSUM.Patient_ID == i].CUSUM.plot(marker = '*', ms = 10, color = color[8], lw = 3), 
    pp = list(Patient_Info[Patient_Info.patient_id == i][['age', 'sex', 'Death']].values[0])    
    plt.title('Age=' + str(pp[0]) + ', Gender=' + pp[1], fontsize = 20)
    plt.xticks(fontsize = 20), plt.yticks(fontsize = 20), 
    if c == 231 or c == 232 or c == 233:
        ax.set_xlabel('')
    else:
        plt.xlabel('LoS in Hours', fontsize = 20)
    if c == 231 or c == 234:
        plt.ylabel('CUSUM', fontsize = 20)
    else:
        ax.set_ylabel('')
    c = c + 1
    x = ax.get_xticks()
    y = ax.get_yticks()
    plt.text(x[0], y[-1], fig[f], color = 'blue', fontsize = 20)
    f = f + 1
plt.tight_layout()
plt.savefig('D:\\Simulation Model\\Change Point Analysis\\Data\\CUSUM for Patients with change.png', format='png', dpi=300, bbox_inches="tight")




#CUSUM for Creatinin only
Creatinin = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Creatinin_Only.csv")
columns_C = ['Date','Patient_ID',  'CUSUM']
CU_SUM_Mid_C = pd.DataFrame(columns = columns_C)
Mean_Mid = []
CU_SUM = []
t = 0
for c in Creatinin.Patient_ID.unique():
    Data = Creatinin[Creatinin.Patient_ID == c].dropna()
    Data.index = Data.Date
    mean = Data.Creatinin.dropna().mean()
    Mean_Mid.append(mean)
    CU_SUM.append([0])
    res = pd.DataFrame([Data.index[0], c, 0]).T
    res.columns = columns_C
    CU_SUM_Mid_C = CU_SUM_Mid_C.append(res)
    d = 0
    for i in Data['Creatinin'].values:
        s = CU_SUM[t][-1] + (i - mean)
        CU_SUM[t].append(s)
        res = pd.DataFrame([Data.index[d], c, s]).T
        res.columns = columns_C
        CU_SUM_Mid_C = CU_SUM_Mid_C.append(res)
        d = d + 1
    t = t + 1

CU_SUM_Mid_C.to_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Creatinin_CUSUM.csv")

#Troponin and Creatinin
Troponin_Creatnin = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Troponin_CreatninFinal.csv")
Case = Troponin_Creatnin[Troponin_Creatnin.Patient_ID == 439505]
plt.scatter(x = np.log(Troponin_Creatnin.Troponin), y = np.log(Troponin_Creatnin.Creatinin))
markers = ['.','o','v','^', '<', '>', '1', '2','3','4','8','s','p','*','h','H','+','x','D','d','|','_']
color_list = []
for name, hex in matplotlib.colors.cnames.items():
    color_list.append(name)
plt.figure(figsize=(14, 8))
plt.xticks(fontsize=14)
plt.yticks(rotation = 'vertical', fontsize=14)
plt.xlabel('Time', fontsize=14)
plt.ylabel('CUMSUM', fontsize=14)



#Patient Information
Patient_Info = pd.read_csv("D:\\Simulation Model\\Change Point Analysis\\Data\\Patient Information2.csv")
ACS_ID = Patient_Info[Patient_Info.patient_id.isin(CU_SUM_Mid_T.Patient_ID)]
ACS_ID['epizod_start_date'] = pd.to_datetime(ACS_ID['epizod_start_date'])
ACS_ID = ACS_ID[['patient_id', 'sex', 'condition', 'result']]
ACS_ID_red = ACS_ID.drop_duplicates()
ACS_ID_red = ACS_ID_red[ACS_ID_red.patient_id.isin(CU_SUM_Mid_C.Patient_ID)]
ACS_ID_red = ACS_ID_red.dropna()