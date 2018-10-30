from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='hydmod',
      version='0.0.1',
      description='Hydrological modeling components and models',
      long_description=readme(),
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering :: GIS',
      ],
      keywords='GIS, raster, vector',
      url='https://github.com/konradhafen/hydmod',
      author='Konrad Hafen',
      author_email='konrad.hafen@gmail.com',
      license='GPLv3',
      packages=['hydmod'],
      install_requires=[
          'gdal',
          'numpy',
          'richdem',
          'datetime',
          'pandas',
      ],
      include_package_data=True,
      zip_safe=False)