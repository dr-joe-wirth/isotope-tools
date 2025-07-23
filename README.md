# Isotope tools
## installation
install the software with these commands:

```bash
git clone https://github.com/dr-joe-wirth/isotope-tools
conda env create -f ./isotope-tools/environment.yml
```

then check that it was installed:

```bash
conda activate isotope-tools
./isotope-tools/mass_shift_abundances.py --help
```

## usage
```text
mass_shift_abundances.py [-faonh]

required arguments:
    -f, --formula           [str] empirical chemical formula
    -a, --abundance-table   [file] tab-delimited table with four columns:
                                     * element
                                     * isotope
                                     * mass shift
                                     * abundance
    -o, --out               [file] output filename

optional arguments:
    -n, --num-cpus          [int] the number of cpus for parallel processing (default 1)
    -h, --help              print this message and exit
```

## example
There is an [example abundance table](./example_table.tsv). It contains exactly four tab-delimited columns:
  * element
  * isotope
  * mass shift
  * abundance

If we wanted to get the abundances of specific mass shifts of the compound dimethylsulfoniopropionate, which has the empirical chemical formula of `C5H10O2S`, we could use the following command:

```bash
conda activate isotope-tools
./isotope-tools/mass_shift_abundances.py --formula C5H10O2S --abundance-table ./isotope-tools/example_table.tsv --out ./dmsp_abundances.tsv
```

which would produce an output like this:

```text
mass	shift	abundance
134	0	0.8952601282678857
135	1	0.05746311533467398
136	2	0.04528746387040204
138	4	0.0003208097273203428
137	3	0.002512280504477125
139	5	1.5452914343389072e-05
140	6	8.980872416767209e-07
141	7	3.503661251957822e-08
142	8	1.12613139426148e-09
143	9	3.035830270648303e-11
144	10	5.478821417486263e-13
145	11	5.8839419675213844e-15
146	12	3.514995623865684e-17
147	13	1.0338380215741874e-19
148	14	1.1304857344578178e-22
149	15	6.615249430993873e-26
150	16	2.415438144358673e-29
151	17	5.9171416676064596e-33
152	18	1.0054042986922073e-36
153	19	1.1942803780523092e-40
154	20	9.769034839130106e-45
155	21	5.258987953406369e-49
156	22	1.6810469446243642e-53
157	23	2.4216235581200404e-58
```

where the first column is the atomic mass, the second column is the mass shift from the most stable isotopes, and the third column is the abundance of any species that possess the given mass.

Some chemical formulas are quite large and require evaluating very large numbers of isotopomers. To address this, multithreading is available with the argument: `--num-cpus`.
