#!/usr/bin/env python

from setuptools import setup
from setuptools.command.install import install
from distutils.command.build import build
from subprocess import call
import glob,os

BASEPATH = os.path.dirname(os.path.abspath(__file__))
STRAINCALL_PATH = os.path.join(BASEPATH,'StrainCall')
STRAINCALL_CPP = glob.glob(os.path.join(STRAINCALL_PATH,'*.cpp'))

class StrainCallBuild(build):
    def run(self):
        # run the original build code
        build.run(self)
        # build StrainCall
        build_path = os.path.abspath(self.build_temp)
        cmd = ['g++','-std=c++11','-o',os.path.join(STRAINCALL_PATH,'StrainCall')]
        cmd.extend(STRAINCALL_CPP)

        def compile():
            call(cmd)

        self.execute(compile, [], 'Compiling StrainCall')

        # copy resulting tool to library build folder
        self.mkpath(self.build_lib)

        if not self.dry_run:
            self.copy_file(os.path.join(STRAINCALL_PATH,'StrainCall'),self.build_lib)


class StrainCallInstall(install):
    def initialize_options(self):
        install.initialize_options(self)
        self.build_scripts = None

    def finalize_options(self):
        install.finalize_options(self)
        self.set_undefined_options('build', ('build_scripts', 'build_scripts'))

    def run(self):
        # run original install code
        install.run(self)

        # install StrainCall executables
        self.copy_tree(self.build_lib, os.path.join(self.exec_prefix,'bin'))
            

install_requires = ['scipy >= 0.15.1','numpy >= 1.9.2','ete2 >= 2.3.1','mpi4py >= 1.3.1']
scripts = glob.glob('scripts/*.py')

setup(
    name = 'rbra',
    description = 'a pipeline to assemble strain-level full-length 16S rRNA genes',
    author = 'Feng Zeng',
    author_email = 'zengfeng@xmu.edu.cn',
    version = '0.2.0',
    url = 'https://github.com/homopolymer/RBRA',
    scripts = scripts,
    install_requires = install_requires,   
    cmdclass={
        'build': StrainCallBuild,
        'install': StrainCallInstall,
    }
)
