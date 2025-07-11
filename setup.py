from setuptools import setup, find_packages

setup(
    name='ganflu',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'toml',
        'biopython',
        'samotools',
        'miniprot',
    ],
    include_package_data=True,
    python_requires='>=3.10.0',
    author='Satoshi Kawato',
    author_email='kawato@kaiyodai.ac.jp',
    description='Genome ANnotation for inFLUenza viruses (GANFLU)',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/satoshikawato/ganflu/',
    entry_points={'console_scripts': ['ganflu = ganflu.ganflu:main',] },
    # Add other relevant information
)
