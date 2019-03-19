import setuptools

setuptools.setup(
    name='molsym',
    version='0.1',
    description='Analytic point group algebra for molecular symmetry operations.',
    author='Ole Hüter',
    packages=['molsym'],
    install_requires=['bidict'],
    keywords=[
        'molecule', 'symmetry', 'operation', 'point group', 'group theory',
        'representation', 'irreducible', 'reducible', 'character', 'Mulliken',
    ],
    url='https://github.com/oh-fv/molsym',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry'],
)
