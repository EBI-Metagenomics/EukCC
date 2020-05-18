=============================
Frequently asked questions:
=============================

Some issues with EukCC came up in the past and will be documented here.

EukCC fails to launch
============================

Can't run EukCC I get the error: "ImportError: cannot import name 'TreeStyle'
------------------------------------------------------------------------------------

This is caused by the PyQT5 libary, which is dependend on having a grafic driver.
Make first sure you have the PyQt5 library installed and if the error persists
you can try to install `libgl1-mesa-glx` on your operating system or in your
container.

Cwd.c: loadable library and perl binaries are mismatched 
---------------------------------------------------------------------------------
We noticed this warning when using perl from conda. This is caused by compilations of 
perl modules from an older or newer perl library.
it is caused by modules needed for GeneMark-ES. Try debugging by running `gmes_petap.pl` 
in the command line. It should output a help page.

EukCC fails to run
=====================


GeneMark-ES failed on this bin
----------------------------------------------------------------

If GeneMark-ES failed on a MAG you can look into the files 
`workfiles/gmes/gmes.log` and `workfiles/gmes/runGMES.log`. 
They capture the very verbose output of GeneMark-ES.

License missing
#########################
It might be that the license of GeneMark-ES is old or not found. Check that
you have the license in the file `.gm_key` in your home directory.


GeneMark-ES failing to train
###################################
This means that GeneMark-ES failed to complete for this particular bin.
This can be caused by very short or very fragmented MAGs. Most of the time
this will mean that your MAG is not even 50 percent complete and you can most 
likely dicard it.

You can always try to predict proteins unsing another method and then supply 
them directly to EukCC.

`runGMES.log` shows: warning on: /[...]/gm_et_linux_64/gmhmme3
##################################################################
We noticed this issue with some machines (Debian/Ubuntu). This is caused by the 
64 bit version of `gmhmme3` 
segfaulting. One solution is using the 32 bit version of GeneMark-ES (or namely gmhmme3).
See: https://github.com/Finn-Lab/EukCC/issues/3

The 64 Bit version is running well on RedHat Linux 7.4 with linux kernel 3.10. 
If you have issues with a certain setup, you can upen a ticket at GitHub and 
I would extend this list here.
