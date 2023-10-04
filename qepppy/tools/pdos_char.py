#! python

import os
import re


def pdos_char(fname='', kpt_list=[], bnd_list=[], thr=0):
    if not os.path.isfile(fname):
        raise FileNotFoundError(f"File fname:'{fname}' does not exist.")
    if   isinstance(kpt_list, str):
        kl = list(map(int, kpt_list.split(',')))
    elif isinstance(kpt_list, list):
        kl = kpt_list
    else:
        raise ValueError("'kpt_list' must be either a list or a string(with elements eparated by a comma")
    if   isinstance(bnd_list, str):
        bl = list(map(int, bnd_list.split(',')))
    elif isinstance(bnd_list, list):
        bl = bnd_list
    else:
        raise ValueError("'bnd_list' must be either a list or a string(with elements eparated by a comma")

    state_l = []

    k = -1
    b = -1
    ck = False
    ce = False
    cp = False
    with open(fname) as f:
        for line in f:
            line = line.rstrip()
            if ' state #' in line:
                state_l.append(line.split(':')[1])
            if ' k = ' in line:
                k  += 1
                b   = -1
                ce  = False
                kpt = list(map(float, list(filter(None, line.split(' ')))[2:5]))
                if k+1 in kl:
                    ck = True
                    print(f'KPT (#{k + 1:5d}):\t{kpt}')
                else:
                    ck = False
            if ck and (' e(' in line or ' e =' in line):
                if ' e(' in line:
                    off=4
                else:
                    off=2
                ce = True
                b +=1

                el = float(list(filter(None, line.split(' ')))[off])
                if b+1 in bl:
                    print(f'\tBND (#{b + 1:3d}):\t{el} eV')
                else:
                    ce = False

            if cp and ' |psi|^2' in line:
                cp = False
            if ce and (' psi = ' in line or cp):
                cp = True
                l = list(filter(None, re.split(r' +psi = | +\+|\*\[#|\]\+', line)))
                for s, wf in zip(l[1::2], l[0::2]):
                    wf = float(wf)
                    if wf >= thr:
                        print(f'\t\t{state_l[int(s) - 1]}:\t{wf * 100:7.3f}%')

def main():
    import sys
    argc = len(sys.argv)
    if not 2<=argc<=5:
        print('Incorrect use. Pleas pass arguments:'
            "\n\t'fname',"
            "\n\t'kpt_list\t(comma separated)',"
            "\n\t'bnd_list\t(optional) (comma separated)'"
            "\n\t'thr\t(optional) (threshold for printing percantages)'")
        exit()
    if argc==2:
        pdos_char(sys.argv[1])
    elif argc==3:
        pdos_char(sys.argv[1], sys.argv[2])
    elif argc==4:
        pdos_char(sys.argv[1], sys.argv[2], sys.argv[3])
    elif argc==5:
        pdos_char(sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]))

if __name__ == '__main__':
    main()
