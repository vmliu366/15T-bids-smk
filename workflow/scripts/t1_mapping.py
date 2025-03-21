import pymp2rage
import nibabel as nib
import numpy as np

def main(snakemake):
    # Load input images
    uni_img = nib.load(snakemake.input.uni_den)
    inv1_img = nib.load(snakemake.input.inv1)
    inv2_img = nib.load(snakemake.input.inv2)
    
    # Get data arrays
    uni_data = uni_img.get_fdata()
    inv1_data = inv1_img.get_fdata()
    inv2_data = inv2_img.get_fdata()

    # Calculate T1 map using pymp2rage
    t1map = pymp2rage.calculate_t1(
        uni=uni_data,
        inv1=inv1_data,
        inv2=inv2_data,
        beta=snakemake.params.beta,
        TI1=snakemake.params.ti1,
        TI2=snakemake.params.ti2
    )
    
    # Handle potential NaN/inf values
    t1map = np.nan_to_num(t1map, posinf=0, neginf=0)
    
    # Save output with same header as UNI image
    t1_img = nib.Nifti1Image(t1map, uni_img.affine)
    nib.save(t1_img, snakemake.output.t1map)

if __name__ == "__main__":
    main(snakemake)