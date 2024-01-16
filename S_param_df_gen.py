import pandas as pd
import numpy as np

def SP_df_gen(PATH, num_modes, run_param = None, run_param_name = None, run_param_ind = None):
    
    if (run_param == None):
        df = pd.DataFrame(columns = [])
    else:
        df = pd.DataFrame(columns = [run_param_name])
    for i in range(2*num_modes):
        df.insert(loc = 2*i+(1 if (run_param != None) else 0), column = 'Re(S' + str(i+1) + '1)', value = None)
        df.insert(loc = 2*i+(2 if (run_param != None) else 1), column = 'Im(S' + str(i+1) + '1)', value = None)

    for j in range((len(run_param) if (run_param != None) else 1)):
        dfarr = []
        x = []
        if run_param != None:
            x.append(run_param[j])
        for i in range(2*num_modes):
            try:
                dfarr.append(pd.read_csv(PATH + (str(run_param[j]) if (run_param != None) else "") + (run_param_ind if (run_param_ind != None) else "") + str(i + 1) + '1.csv'))
            except:
                print("file s_param" + ((str(run_param[j]) if (run_param != None) else "") + (run_param_ind if (run_param_ind != None) else "") + str(i + 1) + '1.csv') + "doesnt exist")
                if int(i) < 9:
                    data = [{"Frequency (Hz)": 500000000, 'Re( S' + str(i + 1) + '1 )': 0, 'Im( S' + str(i + 1) + '1 )': 0}]
                else:
                    data = [{"Frequency (Hz)": 500000000, 'Re( S' + str(i + 1) + ',1 )': 0, 'Im( S' + str(i + 1) + ',1 )': 0}]
                #virtualdf = pd.DataFrame([500000000, np.nan, np.nan], columns = ["Frequency (Hz)", 'Re( S' + str(i + 1) + '1 )', 'Im( S' + str(i + 1) + '1 )'])
                virtualdf = pd.DataFrame(data)
                dfarr.append(virtualdf)
                #print(dfarr)
            if int(i) < 9:
                x.append(dfarr[i]['Re( S' + str(i + 1) + '1 )'][0])
                x.append(dfarr[i]['Im( S' + str(i + 1) + '1 )'][0])
            else:
                x.append(dfarr[i]['Re( S' + str(i + 1) + ',1 )'][0])
                x.append(dfarr[i]['Im( S' + str(i + 1) + ',1 )'][0])
        df.loc[j] = x

    S2 = []
    for j in range(int(2*num_modes)):
        S2.append(df['Re(S' + str(int(j + 1)) + '1)'][:]**2 + df['Im(S' + str(int(j+1)) + '1)'][:]**2)

    S2b = 0                         # Sum of backscattering squared S-Paraneters
    S2f = 0                         # Sum of the forward scattering squared S-Parameters
    for i in range(num_modes):
        #print(i)
        S2b = S2b + S2[i]
        S2f = S2f + S2[i + num_modes]

    S2sum = S2b + S2f
    print("Sum of S-parameters squared: ", S2sum)
    #print("Amount not in TE10 mode: ", S2sum - 2*S2[num_modes])

    return df, (S2, S2[num_modes], S2b, S2sum)