#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from setuptools import setup, find_packages


short_descr = "Package implementing data structures and algorithms for the manipulation of cellular complexes."
readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')


# find version number in src/openalea/cellcomplex/version.py
version = {}
with open("src/openalea/cellcomplex/version.py") as fp:
    exec(fp.read(), version)


setup_kwds = dict(
    name='openalea.cellcomplex',
    version=version["__version__"],
    description=short_descr,
    long_description=readme + '\n\n' + history,
    author="Guillaume Cerutti, Frederic Boudon, ",
    author_email="guillaume.cerutti@inria.fr, frederic.boudon@inria.fr, ",
    url='https://github.com/VirtualPlants/cellcomplex',
    license='cecill-c',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        ],
    tests_require=[
        "coverage",
        "mock",
        "nbconvert",
        "nose",
        "sphinx",
        ],
    entry_points={},
    keywords='',
    test_suite='nose.collector',
)
# #}
# change setup_kwds below before the next pkglts tag

setup_kwds['entry_points']['wralea'] = ['mesh = openalea.mesh_oalab_wralea']
setup_kwds['entry_points']['oalab.applet'] = ['oalab.applet/mesh = openalea.cellcomplex.mesh_oalab.plugin.applet']
setup_kwds['entry_points']['oalab.plugin'] = ['oalab/oalab.core = openalea.cellcomplex.mesh_oalab.plugin.mimedata']
setup_kwds['entry_points']['oalab.world'] = ['oalab.world/mesh = openalea.cellcomplex.mesh_oalab.plugin.world']

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
