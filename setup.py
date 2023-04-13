from setuptools import setup

setup(
        name = 'CountESS Plugin: Minimap2 Wrapper',
        version = '0.0.1',
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
            'countess>=0.0.1',
        ],
        license_files = ('LICENSE.txt',),
        classifiers = [
            'Intended Audience :: Science/Research',
            'License :: Public Domain',
            'Operating System :: OS Independent',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
)

