from setuptools import setup, find_packages

setup(name='layouteditor-wrapper',
      version='0.7.1',
      description='A Python wrapper for LayoutScript, the scripting interface of Juspertor LayoutEditor.',
      author='Daniel Flanigan',
      author_email='daniel.isaiah.flanigan@gmail.com',
      url='https://www.github.com/danielflanigan/layouteditor-wrapper',
      packages=find_packages(where='source'),
      package_dir={'': 'source'},
      install_requires = [
            'numpy'
      ]
      )
