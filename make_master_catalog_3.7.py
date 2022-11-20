'''
Extract info. from McD 107 inch's fits header for
making the master obs. catalog

by Shih-Yun Tang Nov. 19, 2022
'''

from astropy.io import fits
import os, sys
from tqdm import tqdm
import pandas as pd
from datetime import datetime, date
import numpy as np

import warnings
warnings.filterwarnings('ignore')
# pd.set_option("display.max_rows", None)

GET_HEADER_INFO = ['OBJECT', 'OBSERVAT', 'OBSERVER', 'EXPTIME', 'DARKTIME', 
                    'IMAGETYP', 'DATE-OBS', 'UT', 'RA', 'DEC', 'EQUINOX', 
                    'HA', 'ZD', 'AIRMASS']

#------------

def checkIFin_107inch():
    """check if the .py is under the 107inch dir
    """
    current_dir = os.getcwd().split('/')
    if current_dir[-2] != '107inch':
        sys.exit(
            'THIS PROGRAM ONLY WORKS WHEN RUN UNDER THE 107inch/McD_masterlog_generator/. DIR!')
        
def cal_s2n(h, order=25):
    
    flux0 = h[1].data[0][0]
    sig0  = h[1].data[0][1]
    s2n = np.nanmedian(flux0[order]/sig0[order])
    return s2n

def output(df, filename):
    
    today = date.today()
    today_str = today.strftime("%b-%d-%Y")
    
    # df.to_csv(f'./{filename}{today_str}.csv', index=False)
    df.to_excel(f'../{filename}{today_str}.xlsx', index=False, engine='xlsxwriter')
    
    

def main_masterlog(GET_HEADER_INFO):
    
    checkIFin_107inch()
    
    _current_dir_files = os.listdir()
    # make sure only observation data dir is left (i.e., with 2011_oct)
    obs_dirs = [i for i in _current_dir_files if (len(i) == 8) & ('_' in i)] 
    # sorting
    obs_dirs.sort(key = lambda date: datetime.strptime(date, '%Y_%b'))

    info_master_box = [] # master box to save all info list
    out_column_names = GET_HEADER_INFO.copy() # column names for the info_master_box
    
    # loop through all yyyy_mm dirs
    for yyyy_mm in tqdm(obs_dirs):
        
        _current_dir_files = os.listdir(f'../{yyyy_mm}')
        # make sure only observation night's dir is left (i.e., n1)
        n_x_dirs = [i for i in _current_dir_files if (len(i) <=3) & ('n' in i)] 

        # loop through all nights dirs
        for n_x in n_x_dirs:
            
            _current_dir_files = os.listdir(f'../{yyyy_mm}/{n_x}')
            # make sure only spectra data (.ech) is left
            spectrum_echs = [i for i in _current_dir_files if ('.ech' in i)] 
            
            # loop through all spectra
            for ech_file in spectrum_echs:
                
                ech_path = f'../{yyyy_mm}/{n_x}/{ech_file}' # each spectrum's dir
                h = fits.open(ech_path) # open the fits file
                
                info = [] # list to save header info
                info.insert(0, ech_path) # save first column as the spectum dir
                
                # extract header info  
                for head_n in GET_HEADER_INFO:
                    try:
                        info.insert(len(info), h[0].header[f'{head_n}'])
                    except KeyError:
                        info.insert(len(info), 'NaN')   
                        
                # -- add a date-time column for later sorting
                date_time = f"{h[0].header[f'DATE-OBS']} {h[0].header[f'UT']}"
                info.insert(len(info), date_time)
                
                # -- add a s2n column
                if 'obj' in h[0].header[f'IMAGETYP']:
                    try:
                        s2n = cal_s2n(h)
                        info.insert(len(info), s2n )
                        
                    except (TypeError, IndexError):
                        info.insert(len(info), np.nan)
                else:
                    info.insert(len(info), np.nan)
                
                # -- add program type, e.g., rv, ce
                if len(ech_file) > 10:
                    pro_type = ech_file[5:7]
                    info.insert(len(info), pro_type)
                else:
                    info.insert(len(info), 'NaN')
                
                
                # -- if need extra info., add them at below --
                # info.insert(len(info), ??)
                                          
                info_master_box.append(info)

    # transfer list box into pandas df
    out_column_names.insert(0, 'path')
    out_column_names.insert(len(out_column_names), 'date-time')
    out_column_names.insert(len(out_column_names), 's2n')
    out_column_names.insert(len(out_column_names), 'program_type')
    # -- if need extra info., add their column names below --
    
    # transfer list box into df (table)
    df = pd.DataFrame(data=info_master_box, columns=out_column_names)   

    df['s2n'] = np.round(df['s2n'], 1)

    # sort the table base on the obs time
    df["date-time"] = pd.to_datetime(df["date-time"]) 
    df = df.sort_values(by="date-time", ignore_index=True)

    # clean name column
    df_clean = clean_target_name(df)
    # clean alias names
    df_clean_alias = clean_alias(df_clean)

    # save the result master log
    output(df_clean_alias, f'McD107inch_master_log_up2-{obs_dirs[-1]}_CreatedON-')
    # output(df_clean_alias, f'test-{yyyy_mm}')
    
    return df
    

def clean_target_name(df_obj):
    
    df_obj['OBJECT_clean'] = df_obj['OBJECT'].str.replace(' ', '')
    df_obj['OBJECT_clean'] = df_obj['OBJECT_clean'].str.replace('_', '')
    df_obj['OBJECT_clean'] = df_obj['OBJECT_clean'].str.replace('-', '')
    df_obj['OBJECT_clean'] = df_obj['OBJECT_clean'].str.replace(']', '')
    df_obj['OBJECT_clean'] = df_obj['OBJECT_clean'].str.replace('\\', '')
    df_obj['OBJECT_clean'] = df_obj['OBJECT_clean'].str.lower()
    df_obj['OBJECT_clean'] = df_obj['OBJECT_clean'].str.replace('hubblei4', 
                                                                'hubble4')
    return df_obj

def clean_alias(df):

    name_alias = pd.read_excel('./data_4_make_master/alias_names.xlsx') 
    name_alias_clean = clean_target_name(name_alias) 

    hbc_notnan = ~name_alias_clean['HBC'].isnull()
    other1_notnan = ~name_alias_clean['other1'].isnull()

    hbc_from = name_alias_clean['HBC'][hbc_notnan].to_numpy(), 
    hbc_to = name_alias_clean['OBJECT_clean'][hbc_notnan].to_numpy()

    other_from = name_alias_clean['other1'][other1_notnan].to_numpy(), 
    other_to = name_alias_clean['OBJECT_clean'][other1_notnan].to_numpy()

    for i in range(len(hbc_from)):
        df['OBJECT_clean'] = df['OBJECT_clean'].replace(hbc_from[i], hbc_to[i])
    
    for i in range(len(other_from)):
        df['OBJECT_clean'] = df['OBJECT_clean'].replace(other_from[i], other_to[i])
    
    return df

    

def read_RVTargetObservingHistory():
    remain_columns = ['Target', 'HBC', 'C/W', 'N', 'SEP', 'RA', 'Dec', 'Vmag', 
                      'Kmag', 'SpTyp', 'KPNO 4m', 'IRTF']
    
    target_tab = pd.read_excel(
        './data_4_make_master/RVTargetObservingHistory_06Jan2018.xlsx', 
        sheet_name='Summary', header = 1)  
    target_tab = target_tab[remain_columns]
    
    target_tab = target_tab.rename(columns = {'KPNO 4m': 'KPNO 4m-opt',
                                              'IRTF'   : 'IRTF-CSHELL'})
    
    target_tab['KPNO 4m-opt'].replace(0, np.nan, inplace=True)
    target_tab['IRTF-CSHELL'].replace(0, np.nan, inplace=True)
    
    target_tab['Target_c'] = target_tab['Target'].str.replace(' ', '')
    target_tab['Target_c'] = target_tab['Target_c'].str.lower()
    
    return target_tab


def nir_count():
    target_nir_n = pd.read_csv(
        './data_4_make_master/igrins_target_obs_times_all_clean_duplicate.csv')  
    target_nir_n = target_nir_n[['TAR', 'IGRINS_IR_OBS']]
    
    target_nir_n = target_nir_n.rename(columns = {'TAR':'Target_c',
                                              'IGRINS_IR_OBS': 'IGRINS_obs'})
    
    return target_nir_n

def main_count(df):
    # select only target obj.
    df_obj = df[ (df['IMAGETYP'] == 'object') ]
    
    # df_obj = clean_target_name(df_obj)
    
    # get the number of each object's obs times
    df_obj_counts = df_obj.OBJECT_clean.value_counts()
    
    count_all_df = df_obj_counts.to_frame()
    count_all_df.reset_index(inplace=True)
    count_all_df = count_all_df.rename(columns = {'index':'Target_c',
                                          'OBJECT_clean': 'McD107_obs'})
    
    output(count_all_df, 'McD107inch_target_count_CreatedON-')
    
    target_tab = read_RVTargetObservingHistory()
    target_list = target_tab['Target_c'] .to_numpy()
    
    # locate those real observed target also are interested
    tar_loc = [i for i in range(len(df_obj_counts.keys())) if 
               df_obj_counts.keys()[i] in target_list]
     
    count_df = df_obj_counts[tar_loc]
    count_df = count_df.to_frame()
    count_df.reset_index(inplace=True)
    count_df = count_df.rename(columns = {'index':'Target_c',
                                          'OBJECT_clean': 'McD107_obs'})
    
    target_tab = target_tab.join(count_df.set_index('Target_c'), on='Target_c')
    
    target_nir_n = nir_count()
    target_tab = target_tab.join(target_nir_n.set_index('Target_c'), 
                                 on='Target_c')
    
    output(target_tab, 'YSOrv_target_count_CreatedON-')
    
    
            

if __name__ == "__main__":
    
    df = main_masterlog(GET_HEADER_INFO)
    main_count(df)
    
    
    
    
    
    
    