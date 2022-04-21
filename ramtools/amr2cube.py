import os
from academicpython.tools import betterRun

EXE = "/startrek2nb/chongchong/Sam/test_amr2cube/amr2cube_mod"
# EXE = "/startrek2nb/chongchong/Sam/test_amr2cube/amr2cube_mod3"
# EXE = "/startrek/chongchong/softwares/amr2cube/amr2cube_mod"
# EXE = "/startrek/chongchong/softwares/amr2cube/arm2cube_mod"

TYP_DICT = {'den': 1, 'P': 11, 'vx': 2, 'vy': 3, 'vz': 4, 'Bxl': 5,
            'Byl': 6, 'Bzl': 7, 'Bxr': 8, 'Byr': 9, 'Bzr': 10, 'xHII': 12}

def amr2cube(diri, center, width, lma, fo_fmt="amr2cube/{field}.dat",
                 fields=['den', 'xHII'], is_redo=True):
    """ Run amr2cube for a given output with given parameters

    Args:
        diri: (str) input output directory
        diro: (str) output data directory
        center: (str or list_like) the domain center, in boxlen unit. 'c'
            means (.5, .5, .5)
        width: (float) domain width in boxlen unit
        lma: (int) max level

    Returns:

    """
    if isinstance(center, str):
        if center == 'c':
            center = [.5, .5, .5]
    left = [c - width / 2 for c in center]
    right = [c + width / 2 for c in center]
    params = (f"-xmi {left[0]} -xma {right[0]} -ymi {left[1]} -yma {right[1]} "
              f"-zmi {left[2]} -zma {right[2]}")
    if isinstance(fields, str):
        if fields == 'all':
            fields = TYP_DICT.keys()
    if not os.path.isdir(os.path.dirname(fo_fmt)):
        os.makedirs(os.path.dirname(fo_fmt), exist_ok=1)
    for field in fields:
        typ = TYP_DICT[field]
        # denname = f"{out_dir}/Job{jobid}/out{i:02d}_{field}_l{lma}{suffix}.dat"
        dat = fo_fmt.format(field=field)
        if (not is_redo) and os.path.isfile(dat):
            print(f"{dat} exists. Skipping")
        else:
            cmd = (f"{EXE} -inp {diri}"
                   f" -out {dat} -typ {typ} -lma {lma} {params}")
            betterRun(cmd, prt=True, check=True)
            print(f"{dat} created")
    return 0

