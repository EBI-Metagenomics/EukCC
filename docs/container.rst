Container
=================

Using EukCC with a container is maybe the easiest way of 
getting EukCC to run on any server or machine.

We provide explanations on how to use Docker or Singularity. 

Thec containers contain all dependencies but **GeneMark-ES**!
To install GeneMark-ES you will need to either install ist manually in the container or
follow these steps to put GeneMark-ES into the expected location, which is hard coded by EukCC.


Singularity
------------------
*This  tutorial was tested with singularity 2.6 but should work for later versions as well.*

First we can build the container, this will download and build all
dependencies **except GeneMark-ES**.


.. code-block:: shell

    sudo singularity build  # TODO add path

