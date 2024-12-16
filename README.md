
# Automated Cell Culture Splitter
This repository accompanies the 2024 article "Open-source cell culture automation system with integrated cell counting for passaging microplate cultures". It contains the code and original data used to produce the quantitative figures for the article, as referenced in the data availability statement, and also serves as a hub for links to other resources for reproducing the Automated Cell Culture Splitter.

## What is this?
The Automated Cell Culture Splitter (ACCS) is a system that does automated splitting of cells in 96-well plates using an Opentrons OT-2 pipetting robot and a custom cell counting instrument. The design and software are open-source. Check out the preprint for more context.

## Links
### Manuscript
- [Preprint](TODO)**\[TODO\]** on bioRxiv

### System documentation

- [Manuals](https://drive.google.com/drive/folders/1MgdO0HoPbRsYp-P-zQJTYHPDtVA4qac9?usp=drive_link) for building ACCS and the Cell Counting Imager (Google Drive)
- [CAD models](https://cad.onshape.com/documents/f7aeac8a6b627cec67d3facc) for the Cell Counting Imager and CCI flow cell (OnShape)

### Software repos
- [accs-protocol-software-pub](https://github.com/czbiohub-sf/accs-protocol-software-pub):  
  Framework and tooling for generating the protocol scripts that run on the robot
- [accs-cell-counting-imager-pub](https://github.com/czbiohub-sf/accs-cell-counting-imager-pub):  
  Software that powers the Cell Counting Imager

## Contents of this repository

- Directory `accs23_pub_figs/`:  
  Source code and raw data to reproduce the manuscript figures. Refer to `accs23_pub_figs/README.md` for details.

## Maintainers

The current maintainer of this repository is Greg Courville ([:email:](mailto:greg.courville@czbiohub.org)) of the Bioengineering platform at Chan Zuckerberg Biohub San Francisco.
