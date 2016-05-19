import numpy as np

#data_path = '/Users/k/study/absorber_analysis/work/analysis_IM/absorber_doppler/data/1hr/9maps_1hr_AB2/'

data_path = '/home/p/pen/wkhu/scratch/map_result_gbt_absorb/maps/9maps_1hr_AB_random2/'

sec = ''
#sec = '_secA'
#sec = '_secB'

parts = ['part1','part2','part3','part4','part5','part6','part7','part8','part9']

'''
combine data with 9 parts with doppler correction
'''

for part in parts:
    data_bin1 = np.load(data_path+'/gbt_1hr_bin1'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin1'+sec+'_pointcorr_absorb_'+part+'of9_clean_map_I_800.npy')
    data_bin2 = np.load(data_path+'/gbt_1hr_bin2'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin2'+sec+'_pointcorr_absorb_'+part+'of9_clean_map_I_800.npy')
    data_bin3 = np.load(data_path+'/gbt_1hr_bin3'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin3'+sec+'_pointcorr_absorb_'+part+'of9_clean_map_I_800.npy')
    data_bin4 = np.load(data_path+'/gbt_1hr_bin4'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin4'+sec+'_pointcorr_absorb_'+part+'of9_clean_map_I_800.npy')
    data_bin5 = np.load(data_path+'/gbt_1hr_bin5'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin5'+sec+'_pointcorr_absorb_'+part+'of9_clean_map_I_800.npy')
    
    noise_inv_bin1 = np.load(data_path+'/gbt_1hr_bin1'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin1'+sec+'_pointcorr_absorb_'+part+'of9_noise_inv_diag_I_800.npy')
    noise_inv_bin2 = np.load(data_path+'/gbt_1hr_bin2'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin2'+sec+'_pointcorr_absorb_'+part+'of9_noise_inv_diag_I_800.npy')
    noise_inv_bin3 = np.load(data_path+'/gbt_1hr_bin3'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin3'+sec+'_pointcorr_absorb_'+part+'of9_noise_inv_diag_I_800.npy')
    noise_inv_bin4 = np.load(data_path+'/gbt_1hr_bin4'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin4'+sec+'_pointcorr_absorb_'+part+'of9_noise_inv_diag_I_800.npy')
    noise_inv_bin5 = np.load(data_path+'/gbt_1hr_bin5'+sec+'_pointcorr_absorb_'+part+'of9/fir_1hr_bin5'+sec+'_pointcorr_absorb_'+part+'of9_noise_inv_diag_I_800.npy')
    
    
    data_list = [data_bin1,data_bin2,data_bin3,data_bin4,data_bin5]
    noise_inv_list = [noise_inv_bin1,noise_inv_bin2,noise_inv_bin3,noise_inv_bin4,noise_inv_bin5]
    #define the doppler correction factor in 5 doppler bins 
    freq_correct = [1.00007025917,1.0000321345,0.99998527369,0.99995613753,0.999928132211]
    
    spec_comb = []
    data_combT = []
    noise_inv_combT = []
    #define the frequency bins which is as large as 1/4 of the original bins, the values refer to the central value of each bins
    spec = np.arange(900018310.547,700030517.578,-48828.125/4)
    spec /= 1.e6
    #define the frequency bins' edges' value when will be used in  the np.histogram
    spec_edges = np.arange(900018310.547,700018310.5,-48828.125/4) + 48828.125/8
    spec_edges /= 1.e6
    spec_edges.sort()
    #create the combination of 5 doppler bins' frequency which have been applied doppler correction
    for i in np.arange(5):
        spec_comb.extend(spec/freq_correct[i])
    print np.shape(spec_comb)

    #create the combination of 5 doppler bins' data or noise_inv_diag
    for data in data_list:
        dataT = np.repeat(data.T,4)
        dataT = dataT.reshape((len(data.T),len(data.T[0]),4*len(data.T[0][0])))
        data_combT.append(dataT)
    for noise_inv in noise_inv_list:
        noise_invT = np.repeat(noise_inv.T,4)
        noise_invT = noise_invT.reshape((len(noise_inv.T),len(noise_inv.T[0]),4*len(noise_inv.T[0][0])))
        noise_inv_combT.append(noise_invT)
    
    #initialize the final result array
    sr_resultT = np.zeros(np.shape(data_combT[0]))
    
    print np.shape(data_combT)
    #assign the data or noise from different sky pixels to the frequency which have been shifted according to the doppler correction,use the ∑data_k*noise_inv_i/∑noise_inv_i to weight the data

    for i in np.arange(len(data_combT[0])):
        for j in np.arange(len(data_combT[0][0])):
            wdata_combT = []
            noise_inv = []
            for k in np.arange(5):
                noise_inv.append(noise_inv_combT[k][i][j])
            
            noise_inv = np.concatenate(noise_inv)
            count = np.histogram(spec_comb,spec_edges)
            noise_inv_his = np.histogram(spec_comb,spec_edges,weights = noise_inv)
            noise_inv_mean = noise_inv_his[0]
            for k in np.arange(5):
                wdata_combT.append(data_combT[k][i][j]*noise_inv_combT[k][i][j])
                #noise_combT.append(noise_inv_combT[k][i][j]/mean**2)
            wdata_combT = np.concatenate(wdata_combT)
            n = np.histogram(spec_comb,spec_edges)
            spectra = np.histogram(spec_comb,spec_edges,weights = wdata_combT)
            spectra = spectra[0]/n[0]
            spectra = spectra/noise_inv_mean 
            spectra = spectra.tolist()
            spectra.reverse()
            sr_resultT[i][j] = np.array(spectra)
    
    sr_result = sr_resultT.T
    print np.shape(sr_result)
    np.save("/home/p/pen/wkhu/scratch/map_result_gbt_absorb/maps/9maps_1hr_AB_random2/new/sr_"+part+"_9parts_1hr"+sec+".npy",sr_result)
    
    
