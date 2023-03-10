import argparse
import os
import shutil
from pathlib import Path
import multiprocessing
import re

def split_list(alist, splits=1):
    length = len(alist)
    return [alist[i * length // splits: (i + 1) * length // splits]
            for i in range(splits)]

def getSubjectID(path):
    """
    :param path: path to data file
    :return: return the BIDS-compliant subject ID
    """
    stringList = str(path).split("/")
    indices = [i for i, s in enumerate(stringList) if 'sub-' in s]
    text = stringList[indices[0]]
    try:
        found = re.search(r'sub-(\w+)_', text).group(1)
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

def CopyandCheck(orig, target):
    '''
    This function tries to copy a file with current path "orig" to a target path "target" and 
    checks if the orig file exists and if it was copied successfully to target

    :param orig: full path of original location (with (original) filename in the path)
    :param target: full path of target location (with (new) filename in the path)
    '''
    if os.path.exists(orig):
        shutil.copy(orig, target)
        if not os.path.exists(target):
            raise ValueError(f'failed to copy {orig}')
        else:
            print(f'successfully copied {os.path.basename(orig)} to {os.path.basename(target)} target location')
    else:
        raise Warning(f'file {os.path.basename(orig)} does not exist in original folder!')

def process_samseg(dirs, derivatives_dir, freesurfer_path):

    for dir in dirs:

        ### assemble T1w and FLAIR file lists
        t1w = sorted(list(Path(dir).rglob('*T1w*')))
        flair = sorted(list(Path(dir).rglob('*FLAIR*')))

        t1w = [str(x) for x in t1w]
        flair = [str(x) for x in flair]

        if (len(t1w) != len(flair)) or (len(t1w) <= 1) or (len(flair) <= 1):
            # instead of using assert we use this mechanism due to parallel processing
            # assert len(t1w) == len(flair), 'Mismatch T1w/FLAIR number'
            # we do not check for file corresondance as lists are sorted anyway
            print(f"Fatal Error for {dir}")
            break

        ### perform registartion with both T1w images
        # do the registration in a template folder and distribute its results to BIDS conform output directories later
        # create template folder
        temp_dir = os.path.join(derivatives_dir, f'sub-{getSubjectID(t1w[0])}', 'temp')
        temp_dir_output = os.path.join(temp_dir, "output")
        Path(temp_dir_output).mkdir(parents=True, exist_ok=True)
        # pre-define paths of registered images 
        t1w_reg = [str(Path(x).name).replace("T1w.nii.gz", "space-common_T1w.mgz") for x in t1w]
        flair_reg_field = [str(Path(x).name).replace("FLAIR.nii.gz", "space-common_FLAIR.lta") for x in flair]
        flair_reg = [str(Path(x).name).replace("FLAIR.nii.gz", "space-common_FLAIR.mgz") for x in flair]
        # call SAMSEG 
        os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                    cd {temp_dir}; \
                    mri_robust_template --mov {" ".join(map(str, t1w))} --template mean.mgz --satit --mapmov {" ".join(map(str, t1w_reg))};\
                    ')

        ### co-register flairs to their corresponding registered T1w images
        # initialize an empty list fo timepoint argument for samseg that will be used later
        cmd_arg = []
        # iterate over all timepoints and call SAMSEG
        for i in range(len(flair)):
            # get transformation
            os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                        cd {temp_dir}; \
                        mri_coreg --mov {flair[i]} --ref {t1w_reg[i]} --reg {flair_reg_field[i]};\
                        ')

            # apply transformation
            os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                        cd {temp_dir}; \
                        mri_vol2vol --mov {flair[i]} --reg {flair_reg_field[i]} --o {flair_reg[i]} --targ {t1w_reg[i]};\
                        ')

            # generate timepoint argument for samseg
            cmd_arg.append(f'--timepoint {t1w_reg[i]} {flair_reg[i]}') 

            # generate a derivative folder for each session (BIDS)
            deriv_ses = os.path.join(derivatives_dir, f'sub-{getSubjectID(t1w[i])}', f'ses-{getSessionID(t1w[i])}', 'anat')
            Path(deriv_ses).mkdir(parents=True, exist_ok=True)

        ### run SAMSEG longitudinal segmentation 
        os.system(f'export FREESURFER_HOME={freesurfer_path} ; \
                    cd {temp_dir}; \
                    run_samseg_long {" ".join(map(str, cmd_arg))} --threads 4 --pallidum-separate --lesion --lesion-mask-pattern 0 1 -o output/\
                    ')

        ### copy output files from temp folder to their session folders
        # write paths of output folders of the timepoint (tp) in a list
        tp_folder = sorted(list(str(x) for x in os.listdir(temp_dir_output) if "tp" in str(x)))

        # only continue if more than one timepoint was segmented
        if len(tp_folder) > 1:
            # copy the mean image file
            mean_temp_location = os.path.join(temp_dir, "mean.mgz")
            mean_target_location = os.path.join(derivatives_dir, f'sub-{getSubjectID(t1w_reg[0])}', f'sub-{getSubjectID(t1w_reg[0])}' + '_mean.mgz')
            CopyandCheck(mean_temp_location, mean_target_location)
            # iterate through all the timepoints 
            for i in range(len(tp_folder)):
                # initialize empty lists
                tp_files_temp_path = []
                tp_files_ses_path = []

                # define target path in derivatives (BIDS-conform)
                deriv_ses = os.path.join(derivatives_dir, f'sub-{getSubjectID(t1w_reg[i])}', f'ses-{getSessionID(t1w_reg[i])}', 'anat') 
                
                # write template paths and target paths of files of timepoint i into a list (and prepend subject & session IDs (BIDS convention))
                tp_folder_path = os.path.join(temp_dir_output, tp_folder[i])
                tp_files = os.listdir(tp_folder_path)
                for filename in tp_files:
                    # rename to BIDS
                    tp_files_bids = f'sub-{getSubjectID(t1w_reg[i])}' + '_' + f'ses-{getSessionID(t1w_reg[i])}' + '_' + filename
                    # list template and target paths
                    tp_files_temp_path.append(os.path.join(tp_folder_path, filename))
                    tp_files_ses_path.append(os.path.join(deriv_ses, tp_files_bids))

                # define location of files in template folder
                t1w_reg_temp_location = os.path.join(temp_dir, t1w_reg[i])
                flair_reg_temp_location = os.path.join(temp_dir, flair_reg[i])
                flair_reg_field_temp_location = os.path.join(temp_dir, flair_reg_field[i])

                # define location of files in target folder
                t1w_reg_ses_location = os.path.join(deriv_ses, t1w_reg[i])
                flair_reg_ses_location = os.path.join(deriv_ses, flair_reg[i])
                flair_reg_field_ses_location = os.path.join(deriv_ses, flair_reg_field[i])
                
                # copy files from template output folder to target output folder (one output folder per session)
                # each time there is a check if the file in the temp folder exists and if it  was copied successfully
                # T1w
                CopyandCheck(t1w_reg_temp_location, t1w_reg_ses_location)
                # FLAIR
                CopyandCheck(flair_reg_temp_location, flair_reg_ses_location)
                # FLAIR transformation file
                CopyandCheck(flair_reg_field_temp_location, flair_reg_field_ses_location)
                # copy output files to target folder
                for i in range(len(tp_files_temp_path)):
                    CopyandCheck(tp_files_temp_path[i], tp_files_ses_path[i])
        else:
            print(f'-- Only one or no timepoint in template folder, therefore no longitudinal data copied!')

        # delete the temp folder
        shutil.rmtree(temp_dir)
        if os.path.exists(temp_dir):
            raise ValueError(f'failed to delete the template folder: {temp_dir}')
        else:
            print(f'successfully deleted the template folder: {temp_dir}')
            
# end of function definitions
###########################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Run SAMSEG Longitudinal Pipeline on cohort.')
    parser.add_argument('-i', '--input_directory', help='Folder of derivatives in BIDS database.', required=True)
    parser.add_argument('-n', '--number_of_workers', help='Number of parallel processing cores.', type=int, default=os.cpu_count()-1)
    parser.add_argument('-f', '--freesurfer_path', help='Path to freesurfer binaries.', default='/home/twiltgen/Tun_software/Freesurfer/FS_7.3.2/freesurfer')

    # read the arguments
    args = parser.parse_args()

    # generate derivatives/labels/
    derivatives_dir = os.path.join(args.input_directory, "derivatives/samseg-longitudinal-7.3.2")
    Path(derivatives_dir).mkdir(parents=True, exist_ok=True)
    data_root = Path(os.path.join(args.input_directory))

    dirs = sorted(list(data_root.glob('*')))

    # filter directories to only include sub- folders
    dirs = [str(x) for x in dirs]
    dirs = [x for x in dirs if "sub-" in x]

    files = split_list(dirs, args.number_of_workers)

    # initialize multithreading
    pool = multiprocessing.Pool(processes=args.number_of_workers)
    # creation, initialisation and launch of the different processes
    for x in range(0, args.number_of_workers):
        pool.apply_async(process_samseg, args=(files[x], derivatives_dir, args.freesurfer_path))

    pool.close()
    pool.join()
