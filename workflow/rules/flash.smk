# ==================================================================================
# Input Handling Functions
# ==================================================================================
def get_flash_images(wildcards):
    """Organize files by patterns"""
    nifti_dir = checkpoints.dcm_to_nii.get(subject=wildcards.subject).output[0]
    flash_files = sorted(glob(os.path.join(nifti_dir, "*T2star_Flash3D_ISO*.nii*")))

    # check 
    if not flash_files:
        raise FileNotFoundError(f"No 3D FLASH images found for subject {wildcards.subject}")
    return {
        "flash": flash_files[0]
    }

def get_flash_jsons(wildcards):
    """Get FLASH JSON for a subject"""
    nifti_dir = checkpoints.dcm_to_nii.get(subject=wildcards.subject).output[0]
    json_files = sorted(glob(os.path.join(nifti_dir, "*T2star_Flash3D_ISO*.json")))
    if not json_files:
        raise FileNotFoundError(f"No JSON files found for subject {wildcards.subject}")
    return json_files

# ==================================================================================
# Core Processing Rules
# ==================================================================================

rule flash_reorient_to_ras:
    """Reorient all images to RAS+ space"""
    input:
        flash = lambda wc: get_flash_images(wc)["flash"],
        _nifti_dir=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0]
    output:
        bids_flash = bids(
            root=out_path("bids"),
            subject='{subject}',
            datatype='anat',
            acq='3DFLASH',
            suffix='T1w',
            extension='.nii.gz'
        ),
    shell:
        """
        c3d {input.flash} -orient RIA -o {output.bids_flash} 
        fslswapdim {output.bids_flash}  x z -y {output.bids_flash} 
        """

rule get_bruker_flash_jsons:
    input:
        nii_folder=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        dcm_folder=lambda wc: checkpoints.extract_tar.get(subject=wc.subject).output[0],
        _nifti_dir=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        _dicom_dir=lambda wc: checkpoints.extract_tar.get(subject=wc.subject).output[0]
    params:
        jsons = get_flash_jsons
    output:
        method_json=bids(
            root=out_path('bids'),
            subject='{subject}',
            datatype='anat',
            acq='3DFLASH',
            suffix='T1w',
            extension='.json'
        )
    script:
        '../scripts/extract_bruker_info_mtsat.py'