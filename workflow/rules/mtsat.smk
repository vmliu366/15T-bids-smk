
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

    # re-orient 
    cmds.append(f'c3d {t1w_nii} -orient RIA -o {output.t1w_nii}')
    cmds.append(f'c3d $mtw -orient RIA -o {output.mtw_nii}')
    cmds.append(f'c3d $pdw -orient RIA -o {output.pdw_nii}')

    cmds.append(f'fslswapdim {output.t1w_nii} x z -y {output.t1w_nii}')
    cmds.append(f'fslswapdim {output.mtw_nii} x z -y {output.mtw_nii}')
    cmds.append(f'fslswapdim {output.pdw_nii} x z -y {output.pdw_nii}')


    cmds.append(f'cp {t1w_prefix}.json {output.t1w_json}')
    cmds.append(f'cp {mtpdw_prefix}.json {output.pdw_json}')
    cmds.append(f'cp {mtpdw_prefix}.json {output.mtw_json}')
    return ' && '.join(cmds)






rule import_mtsat_files:
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
        "{params.cmd}"
 

       

