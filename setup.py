from distutils.core import setup

setup(
  name = 'ElMTreeIndex',        
  packages = ['ElMTree'],  
  version = '0.1.4',      
  license='GPL3',       
  description = 'A class for performing knn similarity searches on datasets of chemical formula using the ElMD metric',  
  author = 'Cameron Hagreaves',              
  author_email = 'cameron.h@rgreaves.me.uk', 
  url = 'https://github.com/lrcfmd/ElMTree/',   
  download_url = 'https://github.com/lrcfmd/ElMTree/archive/0.1.0.tar.gz',    
  keywords = ['ChemInformatics', 'Materials Science', 'Machine Learning', 'Materials Representation'],   
  install_requires=[            
          'cython',
          'numpy',
          'tqdm',
          'scipy',
          'ElMD>=0.5.12',
      ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',  
    'Intended Audience :: Developers',     
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3) ',   
    'Programming Language :: Python :: 3',      
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10'
  ],
)