#!/usr/bin/env python

from typing import Union
import multiprocessing as mp
from itertools import product
import csv, getopt, os, re, sys


def _parseArgs() -> Union[tuple[str,str,int],bool]:
    """parses command line arguments
    
    Raises:
        BaseException: missing one or more required arguments
    
    Returns:
        Union[tuple[str,str,int],bool]: (formula,filename,cpus) or False if help requested
    """
    # flags
    FORMULA_FLAGS = ('-f', '--formula')
    TABLE_FLAGS = ('-a', '--abundance-table')
    OUT_FLAGS = ('-o', '--out')
    CPU_FLAGS = ('-n', '--num-cpus')
    HELP_FLAGS = ('-h', '--help')
    
    # opts for getopt
    SHORT_OPTS = FORMULA_FLAGS[0][-1] + ':' + \
                 TABLE_FLAGS[0][-1] + ':' + \
                 OUT_FLAGS[0][-1] + ':' + \
                 CPU_FLAGS[0][-1] + ':' + \
                 HELP_FLAGS[0][-1]
    LONG_OPTS = (FORMULA_FLAGS[1][2:] + '=',
                 TABLE_FLAGS[1][2:] + '=',
                 OUT_FLAGS[1][2:] + '=',
                 CPU_FLAGS[1][2:] + '=',
                 HELP_FLAGS[1][2:])
    
    # helper function to print help message
    def printHelp():
        EOL = '\n'
        GAP = ' '*4
        WIDTH = max([len(x) + len(GAP) + 4 for x in LONG_OPTS])
        
        def getLeftCol(args:tuple[str,str]) -> str:
            return f'{GAP}{", ".join(args):<{WIDTH}}'
        
        print(f'{EOL}gets mass shift abundances for a chemical compound{EOL}' + \
              f'{GAP}Joseph S. Wirth, 2025{EOL*2}' + \
              f'usage:{EOL}' + \
              f'{GAP}{os.path.basename(__file__)} [-{SHORT_OPTS.replace(':', '')}]{EOL*2}' + \
              f'required arguments:{EOL}' + \
              f'{getLeftCol(FORMULA_FLAGS)}[str] empirical chemical formula{EOL}' + \
              f'{getLeftCol(TABLE_FLAGS)}[file] tab-delimited table with four columns:\n'+ \
              f'{GAP}{"":<{WIDTH}}{GAP*2} * element\n' + \
              f'{GAP}{"":<{WIDTH}}{GAP*2} * isotope\n' + \
              f'{GAP}{"":<{WIDTH}}{GAP*2} * mass shift\n' + \
              f'{GAP}{"":<{WIDTH}}{GAP*2} * abundance\n' + \
              f'{getLeftCol(OUT_FLAGS)}[file] output filename{EOL*2}' + \
              f'optional arguments:{EOL}' + \
              f'{getLeftCol(CPU_FLAGS)}[int] the number of cpus for parallel processing (default 1){EOL}' + \
              f'{getLeftCol(HELP_FLAGS)}print this message and exit{EOL}')
    
    # print help if requested
    if len(sys.argv) == 1 or HELP_FLAGS[0] in sys.argv or HELP_FLAGS[1] in sys.argv:
        printHelp()
        return False
    
    # initialize default value
    cpus = 1
    abdFn = None
    outFn = None
    formula = None
    
    # parse the args
    opts,args = getopt.getopt(sys.argv[1:], SHORT_OPTS, LONG_OPTS)
    for opt,arg in opts:
        # formula argument
        if opt in FORMULA_FLAGS:
            formula = arg
        
        # table argument
        elif opt in TABLE_FLAGS:
            abdFn = arg
        
        # outfile argument
        elif opt in OUT_FLAGS:
            outFn = arg
        
        # num cpus argument
        elif opt in CPU_FLAGS:
            cpus = int(arg)
    
    if None in (formula,abdFn,outFn):
        raise BaseException('all required arguments must be specified')
    
    return formula,abdFn,outFn,cpus


def _parseFormula(formula:str) -> dict[str,int]:
    """parses a string chemical formula into a dictionary

    Args:
        formula (str): empirical chemical formula

    Returns:
        dict[str,int]: {atom: count}
    """
    PATTERN = re.compile(r'([A-Z][a-z]{0,1})(\d*)')
    return {atom:int(num) if num else 1 for atom,num in PATTERN.findall(formula)}


def _parseAbundanceTable(fn:str) -> dict[str,dict[int,tuple[int,float]]]:
    out = dict()
    with open(fn, 'r') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            element = row['element']
            isotope = int(row['isotope'])
            shift = int(row['mass shift'])
            abund = float(row['abundance'])
            
            out[element] = out.get(element, dict())
            out[element][isotope] = (shift, abund)
    return out


def _checkCongruency(formula:dict[str,int], abundances:dict[str,dict[int,tuple[int,float]]]) -> None:
    """checks that all atoms in the formula were in the abundance table

    Args:
        formula (dict[str,int]): the dictionary produced by _parseFormula
        abundances (dict[str,dict[int,tuple[int,float]]]): the dictionary produced by _parseAbundanceTable

    Raises:
        BaseException: one more atoms in the formula are missing from the abundance table
    """
    # determine if any atoms are missing from the abundance table
    missing = set(formula.keys()).difference(set(abundances.keys()))
    
    # report missing atoms and quit
    if missing != set():
        raise BaseException(f'missing abundances for the following atoms: {", ".join([f"\'{x}\'" for x in missing])}')


def _getOneIsotopomer(atoms:list[tuple[str,int,int,float]]) -> tuple[int,int,float]:
    """calculates values for a single isotopomer

    Args:
        atoms (list[tuple[str,int,int,float]]): a list of atoms: (element, isotope, mass shift, abundance)

    Returns:
        tuple[str,int,int,float]: the values for the isotopomer: (formula, molecular weight, mass shift, abundance)
    """
    # initialize outputs
    mass = 0
    shift = 0
    abundance = 1.0
    
    # for each atom, atomic mass, mass shift, and abundance in the list
    for iso,ms,ab in atoms:
        # update the output values
        mass += iso
        shift += ms
        abundance *= ab
    
    return mass, shift, abundance


def _getMassShiftAbundances(formula:dict[str,int], abdTable:dict[str,dict[int,tuple[int,float]]], cpus:int) -> dict[tuple[int,int],float]:
    """runs _getOneIsotopomer in parallel and uses the results to get the mass shift abundances

    Args:
        formula (dict[str,int]): the dictionary produced by _parseFormula
        cpus (int): the number of cpus for parallel processing

    Returns:
        dict[tuple[int,int],float]: {(molecular weight, mass shift): abundance}
    """
    # generates arguments for _getOneIsotopomer
    def genArgs():
        expanded_formula = list()
        for atom,num in formula.items():
            expanded_formula.extend([[(iso, shift, ab) for iso, (shift,ab) in abdTable[atom].items()]] * num)
        
        for atoms in product(*expanded_formula):
            yield list(atoms)
    
    # run _getOneIsotopomer in parallel
    with mp.Pool(cpus) as pool:
        try:
            isotopomers = pool.map(_getOneIsotopomer, genArgs())
        except KeyboardInterrupt:
            pool.close()
            pool.join()
            raise KeyboardInterrupt()
    
    # initialize output
    out = dict()
    
    # for each isotopomer
    for mass,shift,abundance in isotopomers:
        # combine abundances for all isotopomers at a each mass shift
        out[(mass, shift)] = out.get((mass, shift), 0)
        out[(mass, shift)] += abundance

    return out


def _writeResults(fn:str, massAbundance:dict[tuple[int,int],float]) -> None:
    """writes results to file

    Args:
        fn (str): the filename to write to
        massAbundance (dict[tuple[int,int],float]): the dictionary produced by _getMassShiftAbundances
    """
    # open the file and create a writer
    with open(fn, 'w') as fh:
        writer = csv.DictWriter(fh, fieldnames=('mass', 'shift', 'abundance'), delimiter='\t')
        writer.writeheader()
        
        # write each row to file
        for (mass,shift),abundance in massAbundance.items():
            writer.writerow({'mass': mass, 'shift': shift, 'abundance': abundance})


def _main() -> None:
    """main runner function
    """
    # parse command line
    args = _parseArgs()
    
    # run if help wasn't requested
    if args:
        formula,abdFn,outFn,cpus = args
        formula = _parseFormula(formula)
        abundance = _parseAbundanceTable(abdFn)
        _checkCongruency(formula, abundance)
        masses = _getMassShiftAbundances(formula, abundance, cpus)
        _writeResults(outFn, masses)


# entrypoint
if __name__ == "__main__":
    _main()
