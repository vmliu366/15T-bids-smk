# grabs inversion images (complex, channel data), and performs mp2rage processing
# TODO:
#  - [x] create UNI-DEN image for each channel
#  - [x] perform channel combination 
#  - [ ] perform T1 mapping  (qMRlab? pymp2rage?)


#mp2rage is  complex conjugate of inv1 * inv 2
# (a-ib)*(c+id)
#  ac +iad - ibc +bd
#  ac+bd +i(bc-ad)
# Re(inv1)*Re(inv2) + Im(inv1)*Im(inv2)  + i ( Im(inv1) * Re(inv2)  - Im(inv1)*Im(inv2)  )
# numerator = Re(above) = Re(inv1)*Re(inv2) + Im(inv1)*Im(inv2)
# denominator = a^2 + b^2 + c^2 + d^2

n_channels=config['n_channels']

wildcard_constraints:
    part='real|imag'


#to motion-correct, use the scanner-provided UNI images 

def get_uni_mp2rage_niftis(wildcards,input):
    niftis = sorted(glob( input.nii_folder + f'/*2_cfmmMP2RAGE_ISO100*.nii'))
    return niftis

def get_corr_uni_mp2rage_niftis(wildcards,input,output):
    in_niftis = sorted(glob( input.nii_folder + f'/*2_cfmmMP2RAGE_ISO100*.nii'))
    out_niftis = [ Path(output.nii_folder) / Path(nii).name for nii in in_niftis ] 
    return out_niftis


def get_uni_mp2rage_jsons(wildcards,input):
    jsons = sorted(glob( input.nii_folder + f'/*2_cfmmMP2RAGE_ISO100*.json'))
    return jsons


def get_flo_indices(wildcards,input):
    niftis = sorted(glob( input.nii_folder + f'/*2_cfmmMP2RAGE_ISO100*.nii'))
    return " ".join([f"{i}" for i in range(1, len(niftis))])

def get_flo_imgs(wildcards,input):
    niftis = sorted(glob( input.nii_folder + f'/*2_cfmmMP2RAGE_ISO100*.nii'))
    return " ".join(niftis[1:])

"""
rule fixorient_uni:
    input:
        nii_folder='niftis/sub-{subject}',
    params:
        in_unis=get_uni_mp2rage_niftis,
        out_unis=get_corr_uni_mp2rage_niftis,
        ref_nii=lambda wildcards, input: sorted(glob(input.nii_folder + f'/*cfmmMP2RAGE_ISO100_01.nii'))[0]
    output: 
        nii_folder=directory('niftis_corrorient/sub-{subject}')
    run:
        shell('mkdir -p {output.nii_folder}')
        for in_uni,out_uni in zip(params.in_unis,params.out_unis):
            shell('{c3d} {in_uni} -swapdim LSA -popas SWAPPED {params.ref_nii} -push SWAPPED -copy-transform -o {out_uni}') 
"""

rule moco_uni:
    """ motion-correct the multiple mp2rage scans using the UNI """
    input:
        nii_folder='niftis/sub-{subject}'
    params:
        UNIs=get_uni_mp2rage_niftis,
        flo_indices=get_flo_indices,
        flo_imgs=get_flo_imgs,
    output:
        affine_dir=directory(
            bids(
                root='work',
                suffix="transforms",
                desc="moco",
                datatype='anat',
                subject='{subject}'
            )
        ),
        nii_4d=bids(
            root='work',
            suffix="UNIs.nii.gz",
            desc="moco",
            datatype='anat',
            subject='{subject}'
        ),
        nii_avg3d=bids(
            root='work',
            suffix="UNI.nii.gz",
            desc="moco",
            datatype='anat',
            subject='{subject}'
        ),
    threads: 8 
    resources:
        mem_mb=32000,
    shadow:
        "minimal"
    group:
        "subj"
    shell:
        """
        dedent () {{
            xargs -L1 echo
        }}
        parallel --eta --jobs {threads} --link \\
            reg_aladin -flo {{2}}  -ref {params.UNIs[0]} -res warped_{{1}}.nii \\
                -aff affine_xfm_ras_{{1}}.txt --rigOnly \\
            ::: {params.flo_indices} \\
            ::: {params.flo_imgs}
             
        mkdir -p {output.affine_dir}
        if [ -z "{params.flo_indices}" ]; then # If no moving volumes, skip motion correction
            echo "Only one volume found. Skipping motion correction."
            # Create dummy transform (identity matrix)
            echo -e '1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 0 1' > {output.affine_dir}/affine_xfm_ras_000.txt
            # Copy the single volume as the 4D and 3D outputs
            cp {params.UNIs[0]} {output.nii_4d}
            cp {params.UNIs[0]} {output.nii_avg3d}
        else
            # Run motion correction for multiple volumes
            parallel --eta --jobs {threads} --link \\
                reg_aladin -flo {{2}} -ref {params.UNIs[0]} -res warped_{{1}}.nii \\
                    -aff affine_xfm_ras_{{1}}.txt --rigOnly \\
                ::: {params.flo_indices} \\
                ::: {params.flo_imgs}
            cp affine_xfm_ras_*.txt {output.affine_dir}
            mrcat {params.UNIs[0]} warped_*.nii {output.nii_4d}
            mrmath {output.nii_4d} mean {output.nii_avg3d} -axis 3
        fi
        """ 


def get_avg_mp2rage_cmd(wildcards,input,output):
    niftis = sorted(glob( input.nii_folder + f'/*cfmmMP2RAGE_ISO100*{wildcards.part}*.nii'))
    xfms = sorted(glob( input.affine_dir + f'/affine_xfm_ras_*.txt'))
    cmds = []

    print({f'/*cfmmMP2RAGE_ISO100*{wildcards.part}*.nii'})
    print(input.nii_folder)

    #create upsampled ref
#    cmds.append(f'c4d {niftis[0]} -resample 200x200x200x100% ref_upsampled.nii.gz && ')

    #apply registration
    ref = niftis[0]
#    ref = 'ref_upsampled.nii.gz'


    for i,(flo,xfm) in enumerate(zip(niftis,xfms)):
        cmds.append(f'reg_resample -ref {ref} -flo {flo} -res resampled_{i}.nii.gz -aff {xfm} && ')

    #average images
    cmds.append(f'{c4d} resampled_*.nii.gz')
    for i in range(len(niftis)-1):
        cmds.append('-add')
    cmds.append(f'-scale {1.0/len(niftis)}')
    cmds.append(f'-o {output}')
    return ' '.join(cmds)




rule avg_mp2rage_complex_channels:
    """ average real or imag scans (no moco here yet, anesthetized anyhow..)"""
    input:
        nii_folder='niftis/sub-{subject}',
        affine_dir=bids(
                root='work',
                suffix="transforms",
                desc="moco",
                datatype='anat',
                subject='{subject}'
            )
    params:
        cmd = get_avg_mp2rage_cmd
    shadow: 'minimal'
    output:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='{part}',suffix='MP2RAGEchannels.nii.gz')
    shell:
        '{params.cmd}'

rule split_inversions:
    input:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='{part}',suffix='MP2RAGEchannels.nii.gz')
    params:
        vols_inv1=f'0:{n_channels-1}',
        vols_inv2=f'{n_channels}:-1'
    output:
        inv1=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='{part}',inv='1',suffix='MP2RAGEchannels.nii.gz'),
        inv2=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='{part}',inv='2',suffix='MP2RAGEchannels.nii.gz')
    shell:
        '{c4d} {input} -slice w {params.vols_inv1} -tile w -o {output.inv1} && '
        '{c4d} {input} -slice w {params.vols_inv2} -tile w -o {output.inv2} '

 
       
rule mp2rage_numerator_denominator:
    input:
        re_inv1=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='real',inv='1',suffix='MP2RAGEchannels.nii.gz'),
        re_inv2=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='real',inv='2',suffix='MP2RAGEchannels.nii.gz'),
        im_inv1=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='imag',inv='1',suffix='MP2RAGEchannels.nii.gz'),
        im_inv2=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='imag',inv='2',suffix='MP2RAGEchannels.nii.gz'),
    output:
        numerator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNInumeratorchannels.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNIdenominatorchannels.nii.gz'),
    shell:
        '{c4d} {input.re_inv1} {input.re_inv2} -multiply '
        ' {input.im_inv1} {input.im_inv2} -multiply '
        ' -add -as NUMERATOR -o {output.numerator} '
        ' {input.re_inv1} -dup -multiply '
        ' {input.im_inv1} -dup -multiply '
        ' {input.re_inv2} -dup -multiply '
        ' {input.im_inv2} -dup -multiply '
        ' -add -add -add -as DENOMINATOR -o {output.denominator} '
   
rule mp2rage_add_bias_term: 
    input:
        numerator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNInumeratorchannels.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNIdenominatorchannels.nii.gz'),
    params:
        num_offset = lambda wildcards: '-shift {offset}'.format(offset=-float(wildcards.beta)),
        den_offset = lambda wildcards: '-shift {offset}'.format(offset=2*float(wildcards.beta)),
    output:
        numerator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNInumeratorchannels.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNIdenominatorchannels.nii.gz'),
    shell:
        '{c4d} {input.numerator} {params.num_offset} -o {output.numerator} && '
        '{c4d} {input.denominator} {params.den_offset} -o {output.denominator} '

rule mp2rage_division_uniden:
    input:
        numerator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNInumeratorchannels.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNIdenominatorchannels.nii.gz'),
    output:
        uni=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNIDENchannels.nii.gz'),
    shell: 
        #c4d divide was not working properly, so using fsl:
        #'c4d {input.numerator} {input.denominator} -divide -replace inf 1000 -inf -1000 NaN 0  -o {output.uni}' 
        'fslmaths {input.numerator} -div {input.denominator} -add 0.5 {output.uni}' #need to add 0.5 to ensure the image has minimum at 0

rule mp2rage_division:
    input:
        numerator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNInumeratorchannels.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNIdenominatorchannels.nii.gz'),
    output:
        uni=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNIchannels.nii.gz'),
    shell: 
        'fslmaths {input.numerator} -div {input.denominator}  {output.uni}' 


def get_sos_combine_cmd(wildcards,input,output):
    cmds=[]
    cmds.append(f'fslsplit {input} channels_ -t')
    cmds.append('&&')
    cmds.append(f'{c3d} channels_*.nii.gz ')
    for i in range(n_channels):
        cmds.append(f'-dup -multiply -popas S{i}') #square it, add to stack

    for i in range(n_channels):
        cmds.append(f'-push S{i}') #push back on stack

    #now we have all the squared channels on the stack, add them up
    for i in range(n_channels-1): #adds prev 2 images, so need one less than the total images
        cmds.append('-add') 

    cmds.append(f'-o {output}')
    return ' '.join(cmds)

rule sos_combine_uni:
    input:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNIchannels.nii.gz'),
    params:
        cmd=get_sos_combine_cmd,
    output:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',suffix='UNI.nii.gz'),
    shadow: 'minimal'
    shell:
        '{params.cmd}'




rule sos_combine_uniden:
    input:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNIDENchannels.nii.gz'),
    params:
        cmd=get_sos_combine_cmd,
    output:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta='{beta}',suffix='UNIDEN.nii.gz'),
    shadow: 'minimal'
    shell:
        '{params.cmd}'


rule sos_combine_invs:
    input:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='{part}',inv='{inv}',suffix='MP2RAGEchannels.nii.gz')
    params:
        cmd=get_sos_combine_cmd,
    output:
        bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='{part}',inv='{inv}',suffix='MP2RAGE.nii.gz')
    shadow: 'minimal'
    shell:
        '{params.cmd}'

rule real_imag_to_phase_mag:
    input:
        real=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='real',inv='{inv}',suffix='MP2RAGE.nii.gz'),
        imag=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='imag',inv='{inv}',suffix='MP2RAGE.nii.gz')
    output:
        mag=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='mag',inv='{inv}',suffix='MP2RAGE.nii.gz'),
        phase=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',part='phase',inv='{inv}',suffix='MP2RAGE.nii.gz'),
    shell:
        '{c3d} {input.imag} {input.real} -atan2 -o {output.phase} && '
        '{c3d} {input.imag}  -dup -multiply -popas IMAG_SQUARED {input.real} -dup -multiply -popas REAL_SQUARED -push IMAG_SQUARED -push REAL_SQUARED -add -sqrt -o {output.mag}'


rule get_bruker_params:
    input:
        nii_folder='niftis/sub-{subject}',
        dcm_folder='dicoms/sub-{subject}'
    params:
        jsons = get_uni_mp2rage_jsons
    output:
        method_json=bids(root='work',subject='{subject}',datatype='anat',acq='mp2rage',suffix='MP2RAGE.method.json'),
    script:
        '../scripts/extract_bruker_info_mp2rage.py'

rule fixorient_uni:
    input:
        uni=bids(
            root='work',
            suffix="UNI.nii.gz",
            desc="moco",
            datatype='anat',
            subject='{subject}'
        ),
        ref=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta=config['beta'],suffix='UNIDEN.nii.gz'),
    output: 
        uni=bids(
            root='work',
            suffix="UNI.nii.gz",
            desc="fixorient",
            datatype='anat',
            subject='{subject}'
        ),
    shell:
        '{c3d} {input.uni} -swapdim LSA -popas SWAPPED {input.ref} -push SWAPPED -copy-transform -o {output}'


rule cp_uni_to_bids:
    input:
        uni=bids(root='work',subject='{subject}',acq='mp2rage',datatype='anat',beta=config['beta'],suffix='UNIDEN.nii.gz'),
    output:
        uni=bids(root='bids',subject='{subject}',acq='mp2rage',datatype='anat',suffix='T1w.nii.gz'),
    shell:
        'cp {input} {output}' 
 