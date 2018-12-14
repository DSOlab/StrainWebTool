from setuptools import setup, find_packages

setup(
    name='StrainWebTool',
    version='0.1-b',
    description='WebTool for StrainTool',
    url='https://github.com/DSOlab/StrainWebTool.git',
    author='Xanthos Papanikolaou, Dimitris Anastasiou',
    author_email='xanthos@mail.ntua.gr, danast@mail.ntua.gr',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['Flask>=0.10']
)
