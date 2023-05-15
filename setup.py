from setuptools import setup

setup(
        name = 'CountESS Plugin: Minimap2 Wrapper',
        version = '0.0.3',
        author = 'Nick Moore',
        maintainer = 'Nick Moore',
        maintainer_email = 'nick@zoic.org',
        packages = [ '.' ],
        entry_points = {
            'countess_plugins': [
                'minimap2 = countess_minimap2:MiniMap2Plugin',
            ],
        },
        install_requires = [
            'countess>=0.0.26,<0.1',
            'mappy~=2.24',
        ],
        license_files = ('LICENSE.txt',),
        classifiers = [
            'Intended Audience :: Science/Research',
            'Operating System :: OS Independent',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
)

