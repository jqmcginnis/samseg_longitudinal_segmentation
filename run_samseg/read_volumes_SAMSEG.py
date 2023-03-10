import argparse
import os
from pathlib import Path
import pandas as pd
import re
import numpy as np

###################################################
# define functions

def getSubjectID(path):
    """
    :param path: path to data file
    :return: return the BIDS-compliant subject ID
    """
    stringList = str(path).split("/")
    indices = [i for i, s in enumerate(stringList) if 'sub-' in s]
    text = stringList[indices[0]]
    try:
        found = re.search(r'sub-m(\d{6})', text).group(1)
    except AttributeError:
        found = ''
    return found

def getSessionID(path):
    """
    :param path: path to data file
    :return: return the BIDS-compliant session ID
    """
    stringList = str(path).split("/")
    indices = [i for i, s in enumerate(stringList) if '_ses-' in s]
    text = stringList[indices[0]]
    try:
        found = re.search(r'ses-(\d{8})', text).group(1)
    except AttributeError:
        found = ''
    return found

def getSegList(path):
    '''
    This function lists all "*_seg.mgz"-files that are in the given path. 

    :param path: path to BIDS derviatives database
    :return: return the lists of "*_seg.mgz"-files.
    '''
    seg_ls = sorted(list(Path(path).rglob('*_seg.mgz')))
    return seg_ls

def combineStats(path, subID, sesID):
    '''
    This function reads the _samseg.stats file and the _sbtiv.stats file of a session and merges both into a dataframe

    :param path: path to BIDS derviatives database
    :param subID: the subject ID of the session
    :param sesID: the session ID of the session
    :return: return the volumes in a dataframe
    '''
    # define path of stat files and load them as dataframe
    samseg_path = os.path.join(path, "sub-m"+subID, "ses-"+sesID, "anat", "sub-m"+subID+"_ses-"+sesID+"_samseg.stats")
    tiv_path = os.path.join(path, "sub-m"+subID, "ses-"+sesID, "anat", "sub-m"+subID+"_ses-"+sesID+"_sbtiv.stats")
    df_samseg_stat = pd.read_csv(samseg_path, header=None, names=["ROI", "volume", "unit"])
    df_tiv_stat = pd.read_csv(tiv_path, header=None, names=["ROI", "volume", "unit"])

    # combine _samseg and _sbtiv and clean ROI names
    df = pd.concat([df_samseg_stat, df_tiv_stat])
    df["ROI"]=df["ROI"].str.replace("# Measure ","")
    # transpose dataframe and clean indices and column names
    df = df.loc[:,["ROI", "volume"]].reset_index().drop("index", axis=1)
    df = df.transpose()
    df.columns = list(df.iloc[0,0:])
    df = df.drop(index = 'ROI').reset_index().drop("index", axis=1)
    # add subject ID and session ID to the dataframe
    df["sub-ID"] = "m"+subID
    df["ses-ID"]= sesID
    df_IDs=df[["sub-ID", "ses-ID"]]
    df.drop(labels=["sub-ID", "ses-ID"], axis=1, inplace=True)
    df.insert(0, "ses-ID", df_IDs["ses-ID"])
    df.insert(0, "sub-ID", df_IDs["sub-ID"])

    return df

####################################################
# main script

parser = argparse.ArgumentParser(description='Read Volumes of SAMSEG Longitudinal Segmentation.')
parser.add_argument('-i', '--input_directory', help='Folder of derivatives in BIDS database.', required=True)
parser.add_argument('-o', '--output_directory', help='Destination folder for the output table with volume stats.', required=True)

# read the arguments
args = parser.parse_args()

# define path of the derivatives folder
derivatives_dir = os.path.join(args.input_directory, "derivatives/samseg-longitudinal-7.3.2")

# get a list with the paths of the _seg.mgz files 
# (we do this because we only want subject- and session-IDs for the cases that have been successfully segmented by SAMSEG)
seg_list = getSegList(derivatives_dir)
# initialize empty dataframe into which we will write the stats data of all cases
df_stat = pd.DataFrame()
for i in range(len(seg_list)):
    # get subject and session ID
    subjectID = getSubjectID(seg_list[i])
    sessionID = getSessionID(seg_list[i])
    # get stats of current session
    loop_stat = combineStats(derivatives_dir, subjectID, sessionID)
    # write stats in the final dataframe
    df_stat = pd.concat([df_stat, loop_stat])

# write stats table to .csv file in chosen output directory
df_stat.to_csv(os.path.join(args.output_directory, "volume_stats.csv"), index=False)

################
# convert the volume_stats.csv file to flattened version 
# (e.g., the data of different timepoints of the same subject are next to each other and no more one above/below the other)

# initialize empty dataframe
df_vol1 = pd.DataFrame(columns = df_stat.columns)
df_vol2 = pd.DataFrame(columns = df_stat.columns)
# get all sub-IDs that ar ein the volume_stats file
sub_ls = []
[sub_ls.append(x) for x in list(df_stat["sub-ID"]) if x not in sub_ls]

for i in range(len(sub_ls)):
    loop_subID = sub_ls[i]
    loop_vol = df_stat[df_stat["sub-ID"]==loop_subID]
    # create two separate dataframes for timepoint 1 and timepoint 2
    df_vol1 = pd.concat([df_vol1, loop_vol[loop_vol["ses-ID"]==np.min(loop_vol["ses-ID"])]])
    df_vol2 = pd.concat([df_vol2, loop_vol[loop_vol["ses-ID"]==np.max(loop_vol["ses-ID"])]])

# add label of timepoint 1 or timepoint 2 to column names (keep sub-ID without label)
df_vol1 = df_vol1.add_suffix(".t1").rename(columns={"sub-ID.t1":"sub-ID"})
df_vol2 = df_vol2.add_suffix(".t2").rename(columns={"sub-ID.t2":"sub-ID"})

# merge the two dataframes into one 
# (now the data of different timepoints of the same subject are next to each other and no more one above/below the other)
df_stat_flat = pd.merge(df_vol1, df_vol2, how = 'inner', on = 'sub-ID')

# write the csv file
df_stat_flat.to_csv(os.path.join(args.output_directory, "volume_stats_flat.csv"), index=False)
    

