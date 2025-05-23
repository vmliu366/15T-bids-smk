# grabs inversion images and UNI, performs mp2rage processing
# TODO:
#  - [x] re-orient and split images  
#  - [ ] create UNI-DEN image from UNI
#  - [ ] perform T1 mapping  (qMRlab? pymp2rage?)


# ==================================================================================
# Input Handling Functions
# ==================================================================================
def get_mp2rage_images(wildcards):
    """Organize MP2RAGE files by patterns"""
    nifti_dir = checkpoints.dcm_to_nii.get(subject=wildcards.subject).output[0]
    inv_files = sorted(glob(os.path.join(nifti_dir, "*2_cfmmMP2RAGE_ISO100*.nii*")))
    uni_files = sorted(glob(os.path.join(nifti_dir, "*1_cfmmMP2RAGE_ISO100*.nii*")))

    print(f"Files in directory: {os.listdir(nifti_dir)}")

    # check 
    if not inv_files:
        raise FileNotFoundError(f"No inversion images found for subject {wildcards.subject}")
    if not uni_files:
        raise FileNotFoundError(f"No UNI images found for subject {wildcards.subject}")
    
    return {
        "inv": inv_files[0],  # Using first matching inversion file, should have only 1 anyways
        "uni": uni_files[0]   
    }

def get_split_inversion_cmd(wildcards, input, output):
    """Generate splitting and reorientation commands"""

    print(f"{input}")
    
    return (
        f"fslsplit {input.inv} {wildcards.subject}_vol_ -t && "
        f"mv {wildcards.subject}_vol_0000.nii.gz {output.inv1} && " # rename
        f"mv {wildcards.subject}_vol_0001.nii.gz {output.inv2} "
    )

def get_uni_mp2rage_jsons(wildcards):
    """Get UNI MP2RAGE JSON sidecars for a subject"""
    nifti_dir = checkpoints.dcm_to_nii.get(subject=wildcards.subject).output[0]
    json_files = sorted(glob(os.path.join(nifti_dir, "*2_cfmmMP2RAGE_ISO100*.json")))
    if not json_files:
        raise FileNotFoundError(f"No JSON files found for subject {wildcards.subject}")
    return json_files

# ==================================================================================
# Core Processing Rules
# ==================================================================================
rule split_inversions:
    """Split 4D inversion images into separate INV1/INV2 files"""
    input:
        inv = lambda wc: get_mp2rage_images(wc)["inv"],
        _nifti_dir=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0]
    output:
        inv1 = bids(
            root=out_path("work"),
            subject="{subject}",
            datatype="anat",
            acq="mp2rage",
            inv="1",
            suffix="MP2RAGE",
            extension=".nii.gz"
        ),
        inv2 = bids(
            root=out_path("work"),
            subject="{subject}",
            datatype="anat",
            acq="mp2rage",
            inv="2",
            suffix="MP2RAGE",
            extension=".nii.gz"
        )
    log:
        out_path("logs/split_inversions_{subject}.log")
    params:
        cmd = get_split_inversion_cmd
    shell:
        "{params.cmd}"


rule reorient_to_ras:
    """Reorient all images to RAS+ space"""
    input:
        uni = lambda wc: get_mp2rage_images(wc)["uni"],
        inv1 = rules.split_inversions.output.inv1,
        inv2 = rules.split_inversions.output.inv2
    output:
        bids_inv1 = bids(
            root=out_path("bids"),
            subject="{subject}",
            datatype="anat",
            acq="mp2rage",
            inv="1",
            suffix="MP2RAGE",
            extension=".nii.gz"
        ),
        bids_inv2 = bids(
            root=out_path("bids"),
            subject="{subject}",
            datatype="anat",
            acq="mp2rage",
            inv="2",
            suffix="MP2RAGE",
            extension=".nii.gz"
        ),
        bids_uni = bids(
            root=out_path("bids"),
            subject='{subject}',
            datatype='anat',
            acq='mp2rage',
            desc='UNI',
            suffix='T1w',
            extension='.nii.gz'
        )
    shell:
        """
        c3d {input.uni} -orient RIA -o {output.bids_uni} 
        fslswapdim {output.bids_uni}  x z -y {output.bids_uni}
        c3d {input.inv1} -orient LAS -o {output.bids_inv1}
        fslcpgeom {output.bids_inv1} {output.bids_uni}
        fslswapdim {output.bids_inv1}  -x y z {output.bids_inv1}
        c3d {input.inv2} -orient LAS -o {output.bids_inv2}
        fslcpgeom {output.bids_inv2} {output.bids_uni}
        fslswapdim {output.bids_inv2}  -x y z {output.bids_inv2}
        """

rule get_bruker_params_inversions:
    input:
        nii_folder=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        dcm_folder=lambda wc: checkpoints.extract_tar.get(subject=wc.subject).output[0],
        _nifti_dir=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        _dicom_dir=lambda wc: checkpoints.extract_tar.get(subject=wc.subject).output[0]
    params:
        jsons = get_uni_mp2rage_jsons
    output:
        uni_json=bids(
            root=out_path('bids'),
            subject='{subject}',
            datatype='anat',
            acq='mp2rage',
            desc='UNI',
            suffix='T1w',
            extension='.json'
        ),
        inv1_json=bids(
            root=out_path('bids'),
            subject='{subject}',
            datatype='anat',
            acq='mp2rage',
            inv='1',
            suffix='MP2RAGE',
            extension='.json'
        ),
        inv2_json=bids(
            root=out_path('bids'),
            subject='{subject}',
            datatype='anat',
            acq='mp2rage',
            inv='2',
            suffix='MP2RAGE',
            extension='.json'
        )
    script:
        '../scripts/extract_bruker_info_mp2rage.py'

