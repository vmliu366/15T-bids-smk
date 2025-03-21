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
    inv_files = sorted(glob(os.path.join(nifti_dir, "*1_cfmmMP2RAGE_ISO100*.nii*")))
    uni_files = sorted(glob(os.path.join(nifti_dir, "*2_cfmmMP2RAGE_ISO100*.nii*")))

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
#        f"{c3d} {output.inv1} -orient RAS -o {output.inv1} && " # reorient
#        f"{c3d} {output.inv2} -orient RAS -o {output.inv2}"
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
        {c3d} {input.uni} -orient ALS -o {output.bids_uni} 
        {c3d} {input.inv1} -orient ARS -o {output.bids_inv1}
        fslswapdim {output.bids_inv1} z x y {output.bids_inv1}
        {c3d} {input.inv2} -orient ARS -o {output.bids_inv2}
        fslswapdim {output.bids_inv2} z x y {output.bids_inv2}
        """


rule generate_uniden:
    input:
        uni = rules.reorient_to_ras.output.bids_uni, 
        inv1 = bids(root=out_path('bids'),
                    subject="{subject}",
                    acq="mp2rage",
                    datatype="anat",
                    inv="1",
                    suffix="MP2RAGE",
                    extension=".nii.gz"),
        inv2 = bids(root=out_path('bids'),
                    subject="{subject}",
                    acq="mp2rage",
                    datatype="anat",
                    inv="2",
                    suffix="MP2RAGE",
                    extension=".nii.gz")
    output:
        uni_den = bids(root=out_path('bids'),
                      subject="{subject}",
                      acq="mp2rage",
                      datatype="anat",
                      desc="UNI-DEN",
                      suffix="T1w",
                      extension=".nii.gz")
    params:
        multiplyingFactor=6.0
    script:
        "../scripts/mp2rage_genUniDen_wrapper.py"

rule get_bruker_params:
    input:
        nii_folder=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        dcm_folder=lambda wc: checkpoints.extract_tar.get(subject=wc.subject).output[0],
        _nifti_dir=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        _dicom_dir=lambda wc: checkpoints.extract_tar.get(subject=wc.subject).output[0]
    params:
        jsons = get_uni_mp2rage_jsons
    output:
        method_json=bids(
            root=out_path('bids'),
            subject='{subject}',
            datatype='anat',
            acq='mp2rage',
            desc='UNI',
            suffix='T1w',
            extension='.json'
        )
    script:
        '../scripts/extract_bruker_info_mp2rage.py'


#rule calculate_t1_map:
#   """Calculate T1 map using pymp2rage"""
#    input:
#        uni_den = rules.reorient_to_ras.output.bids_uni_den,
#        inv1 = rules.reorient_to_ras.output.bids_inv1,
#        inv2 = rules.reorient_to_ras.output.bids_inv2
#    output:
#        t1map = bids(
#            root=out_path("bids"),
#            subject="{subject}",
#            datatype="anat",
#            acq="mp2rage",
#            suffix="T1map",
#            extension=".nii.gz"
#        )
#    params:
#        beta = config["beta"],
#        ti1 = config.get("ti1", 0.7),  # Default TI1 (seconds)
#        ti2 = config.get("ti2", 3.5)   # Default TI2 (seconds)
#    conda:
#        "envs/pymp2rage.yaml"  # Recommended for dependency management
#    script:
#        "../scripts/t1_mapping.py"
