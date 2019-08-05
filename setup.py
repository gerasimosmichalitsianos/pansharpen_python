from distutils.core import setup
setup(
  name='pansharpen',
  version='1.0.0',
  scripts=['bin/pansharpen.py',], 
  license='MIT',
  include_package_data=True, 
  long_description=open('README.md').read(),
  entry_points = {
    'console_scripts': ['pansharpen=pansharpen:main'],
  }
)
