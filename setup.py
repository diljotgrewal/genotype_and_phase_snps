import versioneer
from setuptools import setup, find_packages

setup(
    name='genotype_and_phase_snps',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='python utilities for genotyping',
    author='Diljot Grewal',
    author_email='diljot.grewal@gmail.com',
    entry_points={
        'console_scripts': [
            'genotype_and_phase_snps = genotype_and_phase_snps.cli:cli',
        ]
    },
    package_data={'': ['*.py', '*.R', '*.npz', "*.yaml", "data/*", "*.sh"]}
)
