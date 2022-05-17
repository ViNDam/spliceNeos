# spliceNeos
This project identify the neoantigens that are splicing-specific, patient-specific and associates with patients' survival at the transcript level. The list of neoantigens are provided in addition of related molecular information. 

## Requirements
* Unix system with conda and Docker installed

conda install python=3.10
pip install mhcflurry
mhcflurry-downloads fetch models_class1_presentation
pip install biopython
conda install -c bioconda samtools
conda install -c conda-forge r-base
pip install pymp-pypi
pip install scikit-learn scipy pandas matplotlib statannot adjusttext seaborn lifelines
docker pull combinelab/salmon:latest
docker pull fred2/optitype

### funSim
```
git clone https://github.com/thecodingdoc/GOUtil.git
chmod +wx mdsSpectralCombo.py
Add () to mdsSpectralCombo.py line 158 #if using python3

chmod +x mdsSemSim.py
Add () to mdsSemSim.py line 71 #if using python3

The program uses the C++ Boost library, which should be installed on your system before attempting to compile the source code. The Boost library can be found at:

http://www.boost.org
cd path/to/codes/GOUtil
g++ -O3 -o enrich enrich.C utilities.C --std=gnu++11
```

## USAGE
python 
