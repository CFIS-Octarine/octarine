import sys
from daomop.__version__ import version
from setuptools import setup, find_packages

dependencies = ['requests >= 2.7',
                'astropy >= 0.2.5',
                'vos >= 2.0',
                'numpy >= 1.6.1',
                'sip_tpv',
                'Polygon2',
                'scipy',
                'mp_ephem',
                'ginga']


if sys.version_info[0] > 2:
    print('The MOP package is only compatible with Python version 2.7+, not yet with 3.x')
    sys.exit(-1)

console_scripts = ['daomop_populate = daomop.populate:main',
                   'daomop_stationary = daomop.stationary:main',
                   'daomop_cat = daomop.build_cat:main',
                   'daomop_canfar_job = daomop.canfar_job:main',
                   'hpx_map = daomop.hpx_map:main']

setup(name='daomop',
      version=version,
      url='http://github.com/ijiraq/daomop',
      author='''JJ Kavelaars (jjk@uvic.ca),
              Michele Bannister (micheleb@uvic.ca),
              Austin Beauchamp (austinbeauch@gmail.com)''',
      maintainer='JJ Kavelaars',
      maintainer_email='jjk@uvic.ca',
      description="Dominion Astrophysical Observatory Moving Object Pipeline: daomop",
      long_description='See github repp',
      classifiers=['Intended Audience :: Science/Research',
                   'Topic :: Scientific/Engineering :: Astronomy',
                   'Development Status :: 4 - Beta',
                   'Programming Language :: Python :: 2 :: Only',
                   'Operating System :: MacOS :: MacOS X',
                   'Environment :: X11 Applications',
                   'License :: OSI Approved :: GNU General Public License (GPL)',
                   ],
      install_requires=dependencies,
      entry_points={'console_scripts': console_scripts},
      packages=find_packages(exclude=['tests']),
      package_data={'daomop': ['config/*']}
      )
