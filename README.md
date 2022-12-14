# Molecule Synthesizer For Searching Small Molecule that Promotes Cardiac Differentiation of Human Pluripotent Stem Cells.

This program is created for iCeMS Cardiomyocyte Machine Learning Project.

Installation
-

You should have your computer installed anaconda or miniconda, or python.
Anaconda is the package of python that make you possible to create virtual environments, install other libraries.
If you have your computer installed anaconda, the libraries for statistics or machine learning will be installed.
Miniconda is small size package.  The packages will not be installed.
If you are interested in learning more about machine learning, I recommend installing anaconda.
If your computer does not have enough space, install miniconda.
If you install from the official Python installation, you will need to use pip to install the necessary libraries. 
It probably won't be a problem, but since I developed using anaconda, it may cause problems.

### Distributions
- Anaconda: https://www.anaconda.com/products/distribution
- Miniconda: https://docs.conda.io/en/latest/miniconda.html
- Python: https://www.python.org

1. Clone the repository
    ```bash
    git clone https://github.com/Shunichi1100921/SynthesizeNewMoleculeProject.git MoleculeSynthesizer
    ```

2. Install the required packages for using.
    #### Anaconda
    ```bash
    cd MoleculeSynthesizer
    conda create --name mol_env python=3.8 numpy
    conda activate mol_env
    ```

    #### Pip
    ```bash
    cd MoleculeSynthesizer
    pip install numpy
    ```

Usage
---

### Setup
```bash
vi settings.py
```

Molecular information is set in settings.py.
For example, when synthesizing molecules with the following molecule set, the contents of settings.py is as follows.
- benzothiazole: F4
- Amide: F1
- Aryl: F5
- Alcohol1: F11
- Alcohol2: F25
- Modifier: F14
```python:settings.py
# True if you want to create all possible molecules.
# This calculation is very time-consuming.
All_Fragment = False

# Set up a fragment for each part.
# Refer to README.md for the setting method.
benzothiazole = 'F4'
amide = 'F1'
aryl = 'F5'
alcohol1 = 'F11'
alcohol2 = 'F25'
modifier = 'F14'
```
Alcohol and modifier can also be assigned hydrogen. 

In this case, set each value to None.
```python:settings.py
# True if you want to create all possible molecules.
# This calculation is very time-consuming.
All_Fragment = False

# Set up a fragment for each part.
# Refer to README.md for the setting method.
benzothiazole = 'F4'
amide = 'F1'
aryl = 'F5'
alcohol1 = None
alcohol2 = None
modifier = None
```
Multiple fragments may be set for each value.  In this case, all possible combinations of molecules are synthesized.
To set multiple fragments, set them in a list as follows.

For example, the following will create 2 x 2 x 1 x 1 x 2 x 2 = 16 different molecule.
```python:settings.py
# True if you want to create all possible molecules.
# This calculation is very time-consuming.
All_Fragment = False

# Set up a fragment for each part.
# Refer to README.md for the setting method.
benzothiazole = ['F4', 'F15']
amide = ['F1', 'F36']
aryl = 'F5'
alcohol1 = 'F2'
alcohol2 = [None, 'F2']
modifier = [None, 'F12']
```
Which group each fragment belongs to is described in fragment_classification.dat.
A description of each group is provided in fragment_classification.md.
Please refer to them when setting up fragments for synthesizing.

If you want to synthesize molecules using all possible fragments, set All_Fragment True.  
This process is very time-consuming.  It synthesizes nearly 2 million molecules. 
```python:settings.py
# True if you want to create all possible molecules.
# This calculation is very time-consuming.
All_Fragment = True
```

### Run synthesize program
Run the main script.
```bash
python __main__.py
```
