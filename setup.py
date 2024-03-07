from setuptools import setup, find_packages

setup(
    name='outersum',
    version='0.1.0',
    author='Chris Su',
    author_email='suzhiyue8@gmail.com',
    description='A package for calculating and analyzing outer sums from genomic data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/Chrisu8/outersum',  # URL to your package's repository
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
        'scipy',
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.6',
    extras_require={
        'dev': [
            'pytest',
            'check-manifest',
        ],
        'test': [
            'coverage',
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
