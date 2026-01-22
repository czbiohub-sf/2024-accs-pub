# Automated Cell Culture Splitter
This repository accompanies the 2024 article ["Open-source cell culture automation system with integrated cell counting for passaging microplate cultures"](https://doi.org/10.1101/2024.12.27.629034). It contains the code and original data used to produce the quantitative figures for the article, as referenced in the data availability statement, and also serves as a hub for links to other resources for reproducing the Automated Cell Culture Splitter.

## What is this?
The Automated Cell Culture Splitter (ACCS) is a system that does automated splitting of cells in 96-well plates using an Opentrons OT-2 pipetting robot and a custom cell counting instrument. The design and software are open-source. Check out the preprint for more context.

## Links
### Manuscript
- [Preprint](https://www.biorxiv.org/content/10.1101/2024.12.27.629034v2) on bioRxiv

### System documentation

- Manuals for building ACCS and the Cell Counting Imager: (see `accs24-pub-supplements/`)
- [CAD models](https://cad.onshape.com/documents/f7aeac8a6b627cec67d3facc) for the Cell Counting Imager and CCI flow cell (OnShape)

### Software repos
- [accs-protocol-software-pub](https://github.com/czbiohub-sf/accs-protocol-software-pub):  
  Framework and tooling for generating the protocol scripts that run on the robot
- [accs-cell-counting-imager-pub](https://github.com/czbiohub-sf/accs-cell-counting-imager-pub):  
  Software that powers the Cell Counting Imager

## Contents of this repository

- Directory `accs24-pub-figs/`:  
  Source code and raw data to reproduce the manuscript figures. Refer to `accs23_pub_figs/README.md` for details.
- Directory `accs24-pub-supplements/`:  
  Supplements published with the manuscript, in PDF form.

## Maintainers

The current maintainer of this repository is [greg.courville@biohub.org](mailto:greg.courville@biohub.org)
