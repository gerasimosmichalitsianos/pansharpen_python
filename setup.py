from distutils.core import setup
setup(
  name='python-pansharpen',
  version='1.0.0',
  scripts=['bin/pansharpen.py',], 
  license='MIT',
  include_package_data=True, 
  long_description=open('README.txt').read(),
  entry_points = {
    'console_scripts': ['python-pansharpen=pansharpen:main'],
  }
)
