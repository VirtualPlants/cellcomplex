#!/usr/bin/env python
# -*- coding: utf-8 -*-

# {# pkglts, pysetup.kwds
# format setup arguments

from setuptools import setup, find_packages


short_descr = "Package implementing data structures and algorithms for the manipulation of cellular complexes."
readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')


def parse_requirements(fname):
    with open(fname, 'r') as f:
        txt = f.read()

    reqs = []
    for line in txt.splitlines():
        line = line.strip()
        if len(line) > 0 and not line.startswith("#"):
            reqs.append(line)

    return reqs

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
    url='https://github.com/Guillaume Cerutti/cellcomplex',
    license='cecill-c',
    zip_safe=False,

    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=parse_requirements("requirements.txt"),
    tests_require=parse_requirements("dvlpt_requirements.txt"),
    entry_points={},
    keywords='',
    test_suite='nose.collector',
)
# #}
# change setup_kwds below before the next pkglts tag

setup_kwds['entry_points']['wralea'] = ['mesh = openalea.mesh_oalab_wralea']
setup_kwds['entry_points']['oalab.applet'] = ['oalab.applet/mesh = openalea.cellcomplex.mesh_oalab.plugin.applet']

# do not change things below
# {# pkglts, pysetup.call
setup(**setup_kwds)
# #}
