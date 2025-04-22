from distutils.core import setup
setup(name='monochalcogen_polarization',
      version='1.1',
      #package_dir={'becqsdr': 'utils'},
      packages=['MXPol'],
      install_requires = [
        'ase>=3.23',
        'matscipy>=1.1.0',
        'matplotlib>=3.9.0',
        'becqsdr>=1.0',
        'mace-torch>=0.3.12',     
   ]
)
