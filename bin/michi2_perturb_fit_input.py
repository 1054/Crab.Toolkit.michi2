#!/usr/bin/env python
# 

import os, sys, re
import numpy as np
import click
from astropy.table import Table


@click.command()
@click.argument('fit_in', type=click.Path(exists=True))
@click.argument('n_realizations', type=int)
@click.argument('out_base', type=str)
def main(
        fit_in, 
        n_realizations, 
        out_base, 
    ):
    
    table = Table.read(fit_in, format='ascii.commented_header')
    table = table[table.colnames[0:3]]
    colflux = table.colnames[1]
    colferr = table.colnames[2]
    fluxnew = np.array([
        np.random.normal(loc = table[colflux][i], scale = table[colferr][i], size = n_realizations) 
        for i in range(len(table))
    ])
    print(f'fluxnew.shape: {fluxnew.shape}')
    for i in range(n_realizations):
        k = i+1
        out_dir = f'{out_base}_{k}'
        if not os.path.isdir(out_dir):
            os.makedirs(out_dir)
        table[colflux] = fluxnew[:, i]
        out_file = os.path.join(out_dir, fit_in)
        table.write(out_file, format='ascii.fixed_width', delimiter=' ', bookend=True, overwrite=True)
        with open(out_file, 'r+') as fp:
            fp.seek(0)
            fp.write('#')
        print(f'Output to {out_file!r}')
    





if __name__ == '__main__':
    main()

