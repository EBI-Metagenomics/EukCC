=============================
Frequently asked questions:
=============================


**Can't run EukCC I get the error: "ImportError: cannot import name 'TreeStyle'**

This is caused by the PyQT5 libary, which is dependend on having a grafic driver.
Make first sure you have the PyQt5 library installed and if the error persists
you can try to install `libgl1-mesa-glx` on your operating system or in your
container.

GeneMark-ES problems
========================

**EukCC fails with message: GeneMark-ES failed on this bin:** 
This means that GeneMark-ES failed to complete for this particular bin.
This can be caused by very short or very fragmented MAGs. Most of the time
this will mean that your MAG is not even 50 percent complete and you can most 
likely dicard it.

You can always try to predict proteins unsing another method and then supply 
them directly to EukCC.


