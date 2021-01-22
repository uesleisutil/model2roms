import subprocess
from datetime import datetime
import os

def help():
    """
    @compile This is a simple script to call for automatic compiling of all
    fortran files necessary to run the soda2roms package. This is turned on in
    main.py with the compileAll=True

    Call this from command line using python compile.py

    Created by Trond Kristiansen
    https://github.com/trondkr/model2roms
    """


def compileallifort():
    logfile = "compile.log"
    if os.path.exists(logfile): os.remove(logfile)
    log = open(logfile, 'a')
    """Start the processes"""
    print("\n")

    print("Compiling barotropic.f90 to create ==> barotropic.so")
    proc = subprocess.Popen(
        'f2py --verbose --fcompiler=intelem -c -m barotropic barotropic.f90 --f90flags="-no-heap-arrays"',
        shell=True, stdout=subprocess.PIPE, )
    stdout_value = proc.communicate()
    log.writelines(repr(stdout_value))

    print("Compiling interpolation.f90 to create ==> interpolation.so")
    proc = subprocess.Popen(
        'f2py --verbose --fcompiler=intelem -c -m interpolation interpolation.f90 --f90flags="-no-heap-arrays"',
        shell=True, stdout=subprocess.PIPE, )
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print("Compiling fill.f90 to create ==> extrapolate.so")
    proc = subprocess.Popen(
        'f2py --verbose --fcompiler=intelem -c -m extrapolate fill.f90 --f90flags="-no-heap-arrays"',
        shell=True, stdout=subprocess.PIPE, )
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    log.close()

    print("Compilation finished and results written to file => %s" % (logfile))
    print("\n===================================================================")


def compileallgfortran():
    logfile = "compile.log"
    if os.path.exists(logfile): os.remove(logfile)
    log = open(logfile, 'a')
    """Start the processes"""
    print("\n")

    #proc = subprocess.Popen('module swap PrgEnv-pgi PrgEnv-gnu', shell=True, stdout=subprocess.PIPE, )
    #stdout_value = proc.communicate()
    #log.writelines(repr(stdout_value))

    #proc = subprocess.Popen('module unload notur', shell=True, stdout=subprocess.PIPE, )
    #stdout_value = proc.communicate()
    #log.writelines(repr(stdout_value))

    print("Compiling barotropic.f90 to create ==> barotropic.so")
    proc = subprocess.Popen('f2py -c -m barotropic barotropic.f90',
                            shell=True, stdout=subprocess.PIPE, )
    stdout_value = proc.communicate()
    log.writelines(repr(stdout_value))

    print("Compiling interpolation.f90 to create ==> interpolation.so")
    proc = subprocess.Popen('f2py -c -m interpolation interpolation.f90',
                            shell=True, stdout=subprocess.PIPE, )
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    print("Compiling fill.f90 to create ==> extrapolate.so")
    proc = subprocess.Popen('f2py -c -m extrapolate fill.f90',
                            shell=True, stdout=subprocess.PIPE, )
    stdout_value = proc.communicate()[0]
    log.writelines(repr(stdout_value))

    log.close()

    print("Compilation finished and results written to file => %s" % logfile)
    print("\n===================================================================")


def compilefortran(compiler):
    if compiler == "gfortran":
        compileallgfortran()

    if compiler == "ifort":
        compileallifort()


if __name__ == "__main__":

    print("Adding LDFALGS required for Python 3")
    proc = subprocess.Popen('export LDFLAGS="-undefined dynamic_lookup -bundle"',
                            shell=True, stdout=subprocess.PIPE, )

    compilefortran("gfortran")
