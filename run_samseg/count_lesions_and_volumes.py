import nibabel as nib
import numpy as np
import argparse
import os
import sys
from scipy import ndimage as ndi
from skimage.morphology import binary_dilation

parser = argparse.ArgumentParser()

# required
parser.add_argument('-bl', '--baseline', help='Baseline SAMSEG segmentation.', required=True)
parser.add_argument('-fu', '--followup', help='Follow-up SAMSEG segmentation.', required=True)
parser.add_argument('-o', '--output', help='Output directory.', required=True)
# optional
parser.add_argument('--save-images', action='store_true', default=False, help='Save images.')
parser.add_argument('--min-size', type=float, default=15.0, help='Minimum voxel size, in mm^3.')
parser.add_argument('--connectivity', type=int, default=18, help='Connected component connectivity (26 - 18 - 6).')  # 18 as default as in Commowick2018 (MSSeg challenge)
parser.add_argument('--max-overlap', type=float, default=0.3, help='Maximum overlap between a dilated lesion and another existing lesion to classify it as new/disappearing (in percentage of its volume).')
parser.add_argument('--debug', action='store_true', default=False, help='Verbose option, useful for debugging.')

args = parser.parse_args()

# Set-up connectivity matrix
if args.connectivity == 26:
    connectivity = np.ones([3, 3, 3])
elif args.connectivity == 18:
    connectivity = np.ones([3, 3, 3])
    connectivity[0, 0, 0] = 0
    connectivity[-1, -1, -1] = 0
    connectivity[0, -1, -1] = 0
    connectivity[-1, 0, 0] = 0
    connectivity[0, 0, -1] = 0
    connectivity[0, -1, 0] = 0
    connectivity[-1, 0, -1] = 0
    connectivity[-1, -1, 0] = 0
elif args.connectivity == 6:
    connectivity = np.zeros([3, 3, 3])
    connectivity[1, 1, 1] = 1
    connectivity[0, 1, 1] = 1
    connectivity[2, 1, 1] = 1
    connectivity[1, 0, 1] = 1
    connectivity[1, 2, 1] = 1
    connectivity[1, 1, 0] = 1
    connectivity[1, 1, 2] = 1
else:
    print("Invalid value for argument connectivity, exit")
    sys.exit()

# Load baseline
baseline = np.array(nib.load(args.baseline).get_fdata() == 99, dtype=float)  # Assuming Samseg segmentation 
baseline_affine = nib.load(args.baseline).affine
# Compute voxel size
voxelsize_baseline = np.prod(np.sum(baseline_affine[:3, :3] ** 2, axis=0) ** 0.5)
voxel_resolution_baseline = np.sum(baseline_affine[:3, :3] ** 2, axis=0) ** 0.5
# Load follow-up
followup = np.array(nib.load(args.followup).get_fdata() == 99, dtype=float) # Assuming Samseg segmentation
# First count baseline lesions with connected components
lesions_baseline, number_of_lesions_baseline = ndi.label(baseline, connectivity)
# Remove lesions that are smaller than args.min_size
lesion_sizes = np.bincount(lesions_baseline.ravel())[1:] * voxelsize_baseline  # remove background
ok_sizes_baseline = lesion_sizes > args.min_size
# 
thresholded_baseline_lesions = int(np.sum(ok_sizes_baseline))
# Baseline volume
baseline_volume = np.sum(lesion_sizes[ok_sizes_baseline])

# Do the same for followup lesions
lesions_followup, number_of_lesions_followup = ndi.label(followup, connectivity)
# Remove lesions that are smaller than args.min_size
lesion_sizes = np.bincount(lesions_followup.ravel())[1:] * voxelsize_baseline  # remove background
ok_sizes_followup = lesion_sizes > args.min_size
# 
thresholded_followup_lesions = int(np.sum(ok_sizes_followup))
# Followup volume
followup_volume = np.sum(lesion_sizes[ok_sizes_followup])

# Followup - Baseline (i.e., lesion increase)
fu_min_bl = np.array(followup - baseline > 0, dtype=float)  
# First count lesions with connected components
lesions_fu_min_bl, number_of_lesions_fu_min_bl = ndi.label(fu_min_bl, connectivity)
# Remove lesion that are smaller than args.min_size
lesion_sizes = np.bincount(lesions_fu_min_bl.ravel())[1:] * voxelsize_baseline  # remove background
# First discard lesions that are definitevely smaller than args.min_size
# (avoiding removing ring shape differences which might be smaller than args.min_size but filled they are not)
ok_sizes_fu_min_bl = lesion_sizes > 0.7 * args.min_size
# Also compute distance_map (euclidean) to decide if new lesion has an acceptable shape
distance_map = ndi.distance_transform_edt(lesions_fu_min_bl, sampling=voxel_resolution_baseline)
#
enlarging_lesions = 0
new_lesions = 0
for lesion_number in range(1, len(np.unique(lesions_fu_min_bl))):
    
    if not ok_sizes_fu_min_bl[lesion_number - 1]:
        continue

    lesion = lesions_fu_min_bl == lesion_number

    volume = np.sum(lesion) * voxelsize_baseline  # Assuming same resolution for followup image

    if args.debug:
        print("Lesion: " + str(lesion_number))
        print("Volume [mm^3]: " + str(volume))

    # Check if diff lesion has an hole (i.e., it's enlarging)
    # Right now I'm just taking the diff lesion, fill the holes and check if the enlargment is bigger than (TODO: 5%?) than the filled part
    # TODO: not sure if it is really necessary to check for this 5% increase, although ...
    # if we have a mislabeled voxel inside the lesion, we might classify the lesion as enlarging rather than something else
    # Note that I'm also checking if *the filled lesion* fits the min size criteria 
    filled = ndi.binary_fill_holes(lesion)
    filled_volume = np.sum(filled)
    if filled_volume - np.sum(lesion) > 0.05 * filled_volume and filled_volume > args.min_size:
        if args.debug:
            print("Enlarging lesion")
        enlarging_lesions += 1 
        continue

    if volume <= args.min_size:
        if args.debug:   
            print("Remove lesion, too small")
        ok_sizes_fu_min_bl[lesion_number - 1] = False
        continue

    # Finally, check if we have a new solitary lesion or a lesion abutting from another lesion:
    # 1) has roughly a spheroid shape
    # 2) if dilated the overlap with other lesions is small (less than args.max_overlap of its volume)
    dilated_lesion = binary_dilation(lesion)
    overlap = np.sum(np.logical_and(dilated_lesion, baseline) > 0)  
    if args.debug:
        print("Overlap: " + str(overlap))   
    if np.max(distance_map[lesion]) > 1.1 * voxelsize_baseline and overlap < args.max_overlap * volume:
        if args.debug:
            print("New solitary or abutting lesion")
        new_lesions += 1
        continue

    if args.debug:
        print("Remove lesion, it doesn't fit criteria")
    ok_sizes_fu_min_bl[lesion_number - 1] = False

# Lesion increase volume
fu_min_bl_volume = np.sum(lesion_sizes[ok_sizes_fu_min_bl])

# Baseline - Followup (i.e., lesion decrease)
bl_min_fu = np.array(baseline - followup > 0, dtype=float)
# First count lesions with connected components
lesions_bl_min_fu, number_of_lesions_bl_min_fu = ndi.label(bl_min_fu, connectivity)
# Remove lesion that are smaller than args.min_size
lesion_sizes = np.bincount(lesions_bl_min_fu.ravel())[1:] * voxelsize_baseline  # remove background
# First discard lesions that are definitevely smaller than args.min_size
# (avoiding removing ring shape differences which might be smaller than args.min_size but filled they are not)
ok_sizes_bl_min_fu = lesion_sizes > 0.7 * args.min_size
# Also compute distance_map (euclidean) to decide if disappearing lesion has an acceptable shape
distance_map = ndi.distance_transform_edt(lesions_bl_min_fu, sampling=voxel_resolution_baseline)
#
shrinking_lesions = 0
disappearing_lesions = 0
for lesion_number in range(1, len(np.unique(lesions_bl_min_fu))):
    
    if not ok_sizes_bl_min_fu[lesion_number - 1]:
        continue

    lesion = lesions_bl_min_fu == lesion_number

    volume = np.sum(lesion) * voxelsize_baseline  # Assuming same resolution for followup image

    if args.debug:
        print("Lesion: " + str(lesion_number))
        print("Volume [mm^3]: " + str(volume))

    # Check if diff lesion has an hole (i.e., it's shrinking)
    # Right now I'm just taking the diff lesion, fill the holes and check if the enlargment is bigger than (TODO: 5%?) than the filled part
    # TODO: not sure if it is really necessary to check for this 5% increase, although
    # if we have a mislabeled voxel inside the lesion, we might classify the lesion as shrinking rather than something else
    # Note that I'm then checking if *the filled lesion* fits the min size criteria 
    filled = ndi.binary_fill_holes(lesion)
    filled_volume = np.sum(filled)
    if filled_volume - np.sum(lesion) > 0.05 * filled_volume and filled_volume > args.min_size:
        if args.debug:
            print("Shrinking lesion") 
        shrinking_lesions += 1
        continue

    if volume <= args.min_size:
        if args.debug:
            print("Remove lesion, too small")
        ok_sizes_bl_min_fu[lesion_number - 1] = False
        continue

    # Finally, check if we have a disappearing lesion:
    # 1) has roughly a spheroid shape
    # 2) if dilated the overlap with other lesions is small (less than args.max_overlap of its volume)
    dilated_lesion = binary_dilation(lesion)
    overlap = np.sum(np.logical_and(dilated_lesion, followup) > 0) 
    if args.debug:
        print("Overlap: " + str(overlap))   
    if np.max(distance_map[lesion]) > 1.1 * voxelsize_baseline and overlap < args.max_overlap * volume:
        if args.debug:
            print("Disappearing lesion")
        disappearing_lesions += 1
        continue

    # If we are here, we are removing the lesion
    if args.debug:     
        print("Remove lesion, it doesn't fit criteria")
    ok_sizes_bl_min_fu[lesion_number - 1] = False

# Lesion decrease volume
bl_min_fu_volume = np.sum(lesion_sizes[ok_sizes_bl_min_fu])


# Print number of lesions and volumes
print("Number of lesions baseline (bl): " + str(thresholded_baseline_lesions))
print("Total lesion volume: " + str(baseline_volume))

print("Number of lesions followup (fu): " + str(thresholded_followup_lesions))
print("Total lesion volume: " + str(followup_volume))

print("Number of lesions fu - bl: new: " + str(new_lesions) + " enlarging: " + str(enlarging_lesions) + " tot: " + str(enlarging_lesions + new_lesions))
print("Total lesion volume: " + str(fu_min_bl_volume))

print("Number of lesions bl - fu: disappearing: " + str(disappearing_lesions) + " shrinking: " + str(shrinking_lesions) + " tot: " + str(disappearing_lesions + shrinking_lesions))
print("Total lesion volume: " + str(bl_min_fu_volume))

print("Effective number of followup lesions: " + str(thresholded_baseline_lesions + new_lesions - disappearing_lesions))

# Save results to file
np.savez(os.path.join(args.output, "lesion_count_and_total_volume.npz"), bl_les=thresholded_baseline_lesions,
                                                                         bl_vol_les=baseline_volume,
                                                                         fu_les=thresholded_followup_lesions,
                                                                         fu_les_eff=thresholded_baseline_lesions + new_lesions - disappearing_lesions,
                                                                         fu_vol_les=followup_volume,
                                                                         fu_min_bl_les=enlarging_lesions + new_lesions,
                                                                         fu_min_bl_vol_les=fu_min_bl_volume,
                                                                         bl_min_fu_les=disappearing_lesions + shrinking_lesions,
                                                                         bl_min_fu_vol_les=bl_min_fu_volume)

if args.save_images:        

    # Actually mask out small lesions
    # Baseline
    for lesion_number in range(1, len(np.unique(lesions_baseline))):
        if not ok_sizes_baseline[lesion_number - 1]:    
            lesions_baseline[lesions_baseline == lesion_number] = 0
    # Followup
    for lesion_number in range(1, len(np.unique(lesions_followup))):
        if not ok_sizes_followup[lesion_number - 1]:    
            lesions_followup[lesions_followup == lesion_number] = 0   
    # fu - bl
    for lesion_number in range(1, len(np.unique(lesions_fu_min_bl))):
        if not ok_sizes_fu_min_bl[lesion_number - 1]:    
            lesions_fu_min_bl[lesions_fu_min_bl == lesion_number] = 0   
    # bl - fu
    for lesion_number in range(1, len(np.unique(lesions_bl_min_fu))):
        if not ok_sizes_bl_min_fu[lesion_number - 1]:    
            lesions_bl_min_fu[lesions_bl_min_fu == lesion_number] = 0                                    

    # Save images
    img = nib.Nifti1Image(lesions_baseline, baseline_affine)
    nib.save(img, os.path.join(args.output, "bl_lesions.nii.gz"))
    img = nib.Nifti1Image(lesions_followup, baseline_affine)
    nib.save(img, os.path.join(args.output, "fu_lesions.nii.gz"))
    img = nib.Nifti1Image(lesions_fu_min_bl, baseline_affine)
    nib.save(img, os.path.join(args.output, "fu_min_bl_lesions.nii.gz"))
    img = nib.Nifti1Image(lesions_bl_min_fu, baseline_affine)
    nib.save(img, os.path.join(args.output, "bl_min_fu_lesions.nii.gz"))

print("Done!")
