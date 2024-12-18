# accs24-pub-figs
This directory contains the framework and data used to generate the quantitative figures for the ACCS manuscript.

## General info
Code responsible for generating each figure / family of figures is in submodules located in `accs24_pub_figs/figs`.
Output files will be generated in:

 - `$PWD/out_disp`: figure preview images
 - `$PWD/out_print`: full resolution figures
 - `$PWD/out_text`: report output

## External requirements
### Packages
See `environment.yaml`.
### Data
Certain resources such as imaging datasets are not stored in this repo due to file size or license reasons:

 - OpenCell Dragonfly imaging datasets for `PML0435` and `PML0442`: `accs23_pub_figs/figs/zstack_gallery/input_data/PML????/raw_data/*MMStack_*-*-*.ome.tif`.  
 - Font files:
	- `accs23_pub_figs/fonts/HKGrotesk-Regular.ttf`
	- `accs23_pub_figs/fonts/HKGrotesk-Bold.ttf`
	- `accs23_pub_figs/fonts/RobotoMono-Regular.ttf`

Cached z-stacks are included in the repo to allow Fig. 3(C) to generate without downloading several gigabytes of data. The font files in question can be obtained for free from Google here: [Hanken Grotesk](https://fonts.google.com/specimen/Hanken+Grotesk), [Roboto Mono](https://fonts.google.com/specimen/Roboto+Mono)

## Quick start
Set up environment:
```bash
mamba env create -f environment.yml
mamba activate accs24_figgen
```

Generate figures from single submodule:
```bash
python generate.py name_of_module
```
(e.g. `python generate.py cci_norm_fc`)

Generate all figures:
```bash
python generate_all.py
```

## License
The software code in this directory is released under the 3-Clause BSD license -- see BSD-3-Clause.txt

The bundled experimental data is released under the [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International license](https://creativecommons.org/licenses/by-nc-nd/4.0/).
