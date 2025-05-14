"""
Wood, (2018). QUIT: QUantitative Imaging Tools. Journal of Open Source Software, 
3(26), 656, https://doi.org/10.21105/joss.00656
"""
# ==================================================================================
# Input Handling Functions
# ==================================================================================
def get_import_mtsat_cmd(wildcards,input,output):
    t1w_nii = sorted(glob( input.nii_folder + f'/*T1w_Flash3D_ISO100*.nii'))[-1] #get last one only
    mtpdw_nii = sorted(glob( input.nii_folder + f'/*MT_Flash3D_ISO100*.nii'))[-1] #get last one only
    t1w_prefix=t1w_nii[:-4]
    mtpdw_prefix=mtpdw_nii[:-4]

    cmds=[]

    # splitting enhanced DICOM into MT-on (mtw) and MT-off (pdw)
    cmds.append(f'fslsplit {mtpdw_nii} split_mtpd -t')
    cmds.append(f'mtw=`ls split_mtpd* | head -n 1`')
    cmds.append(f'pdw=`ls split_mtpd* | tail -n 1`')

    # re-orient (for some reason MTw/PDw and T1w are in different orientations)
    cmds.append(f'c3d {t1w_nii} -orient RIA -o {output.t1w_nii}')
    cmds.append(f'c3d $mtw -orient LAS -o {output.mtw_nii}')
    cmds.append(f'c3d $pdw -orient LAS -o {output.pdw_nii}')

    cmds.append(f'fslswapdim {output.t1w_nii} x z -y {output.t1w_nii}')
    cmds.append(f'fslcpgeom {output.t1w_nii} {output.mtw_nii}')
    cmds.append(f'fslcpgeom {output.t1w_nii} {output.pdw_nii}')
    cmds.append(f'fslswapdim {output.mtw_nii} -x y z {output.mtw_nii}')
    cmds.append(f'fslswapdim {output.pdw_nii} -x y z {output.pdw_nii}')



    cmds.append(f'cp {t1w_prefix}.json {output.t1w_json}')
    cmds.append(f'cp {mtpdw_prefix}.json {output.pdw_json}')
    cmds.append(f'cp {mtpdw_prefix}.json {output.mtw_json}')
    return ' && '.join(cmds)

def qi_mtr_cmd(wildcards,input,output):
    nii_files = sorted(glob(os.path.join(input.nii_folder, '*MT_Flash3D_ISO100*.nii')))
    if not nii_files:
        raise FileNotFoundError(f"No matching NIfTI files found in {input.nii_folder}")
    mtpdw_nii = nii_files[-1] # select last if there are multiple 
    subject_id = wildcards.subject
    deri_folder = directory(out_path(f'niftis/sub-{subject_id}'))

    cmds=[]
    cmds.append(f"mkdir -p {deri_folder}")
    cmds.append(f"module load quit/3.4")
    cmds.append(f"qi mtr {mtpdw_nii} -v")
    cmds.append(f"mv MTR.nii.gz {output.mtc_nii}")
    cmds.append(f"cp {input.mtc_input_json} {output.mtc_json}")

    return " && ".join(cmds)

# ==================================================================================
# Core Processing Rules
# ==================================================================================
rule split_mtsat_nii:
    input:
        nii_folder=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0]
    params:
        cmd=get_import_mtsat_cmd
    output:
        t1w_nii=bids(root=out_path('bids'),
                subject='{subject}',
                acq='T1w',
                datatype='anat',
                suffix='MTsat.nii.gz'),
        pdw_nii=bids(root=out_path('bids'),
                subject='{subject}',
                acq='PDw',
                datatype='anat',
                suffix='MTsat.nii.gz'),
        mtw_nii=bids(root=out_path('bids'),
                subject='{subject}',
                acq='MTw',
                datatype='anat',
                suffix='MTsat.nii.gz'),
        t1w_json=bids(root=out_path('bids'),
                subject='{subject}',
                acq='T1w',
                datatype='anat',
                suffix='MTsat.json'),
        pdw_json=bids(root=out_path('bids'),
                subject='{subject}',
                acq='PDw',
                datatype='anat',
                suffix='MTsat.json'),
        mtw_json=bids(root=out_path('bids'),
                subject='{subject}',
                acq='MTw',
                datatype='anat',
                suffix='MTsat.json'),
    shadow: 'minimal'
    shell:
        """
        {params.cmd}
        """
 
rule calc_mtc: 
    input:
        nii_folder=lambda wc: checkpoints.dcm_to_nii.get(subject=wc.subject).output[0],
        mtc_input_json=rules.split_mtsat_nii.output.mtw_json
    params:
        cmd=qi_mtr_cmd
    output:
        mtc_nii=bids(root=out_path('derivatives'),
                subject='{subject}',
                acq='3DGRE',
                desc='QUIT',
                suffix='MTRmap',
                extension='.nii.gz'
                ),
        mtc_json=bids(root=out_path('derivatives'),
                subject='{subject}',
                acq='3DGRE',
                desc='QUIT',
                suffix='MTRmap',
                extension='.json'
                ),
    shell:
        """
        {params.cmd}
        """
        
       

