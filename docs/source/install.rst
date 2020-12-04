=============
 Installation
=============

.. _Python: https://www.python.org
.. _NumPy: http://www.numpy.org
.. _SciPy: http://www.scipy.org
.. _Matplotlib: http://matplotlib.org
.. _Jupyter: http://jupyter.org
.. _Jupyter documentation: http://jupyter.readthedocs.io/en/latest/index.html
.. _Pandas: https://pandas.pydata.org
.. _PyWavelets: http://pywavelets.readthedocs.io/en/latest
.. _miniconda: http://conda.io/miniconda.html
.. _anaconda: https://www.anaconda.com/download
.. _conda: http://conda.pydata.org/docs/user-guide/index.html
.. _guide: https://conda.io/docs/user-guide/getting-started.html
.. _attrs: https://www.attrs.org

----------------------------------
For experienced Python/SciPy users
----------------------------------

Requirements
~~~~~~~~~~~~

* Python_ >=3.7
* NumPy_
* SciPy_
* Matplotlib_
* Pandas_
* PyWavelets_
* attrs_ 

Install
~~~~~~~

::

        conda install -c ussl fluxpart

    or:

        pip install fluxpart

Alternatively, to install fluxpart and create a full conda environment, open this link,

`environment.yml <https://raw.githubusercontent.com/usda-ars-ussl/fluxpart/master/conda.recipe/environment.yml>`_

and right click on the page. Select "save as" and name the saved file `environment.yml`.
Navigate to the directory where you saved `environment.yml`, and execute::

    conda env create -f environment.yml
    conda activate fp02


-----------------------------
For novice Python/SciPy users
-----------------------------

**Fluxpart** is a Python 3 module. Using it requires a Python_ language
interpreter as well as a few standard scientific Python libraries.

The recommended way to install **fluxpart** and create the required
computational environment is to use conda_, an open source package and
environment management system that works on Windows, Linux, and Mac OS X.
Conda can be obtained by installing anaconda_, a Python distribution that
includes conda and 100+ useful scientific libraries.
Alternatively, conda can be obtained by installing the smaller miniconda_ 
distribution.  Anaconda can/should be installed at the user level so its
installation and use does not require administrator rights, nor will it affect
other OS-level Python installations already on your computer.
If anaconda (or miniconda) is already installed on your machine, skip the
next paragraph.

Download and install either anaconda_ or miniconda_ for your operating system.
Choose the latest Python 3.x installer (**not** Python 2.7).
Additional instructions can be found on the conda_ installation page.
The conda getting started guide_ is also a useful resource.

To verify that conda is installed, open a terminal window on Linux or OS X or
go to the Start menu on Windows and search for and open
"Anaconda Prompt". At the command prompt, enter::

    conda --version

Conda should respond with the version number of the installed package.
If you have an older, previously installed version of conda,
it is advised to update it  with::

    conda update conda

Next, open this link,

`environment.yml <https://raw.githubusercontent.com/usda-ars-ussl/fluxpart/master/conda.recipe/environment.yml>`_

and right click on the page. Select "save as" and name the saved file `environment.yml`.
Navigate to the directory where you saved `environment.yml`, and execute::

    conda env create -f environment.yml

That will create a conda environment containing **fluxpart** and its dependencies.
Finally, activate the **fluxpart** environment::

    conda activate fp02

The command line prompt should now be prepended with ``(fp02)``,
indicating that the **fluxpart** environment is active in the shell session.

If you are new(ish) to Python_, a good tool for learning and interactively
building-up **fluxpart** analyses is the `Jupyter notebook`__. With the
**fluxpart** environment active, install Jupyter_ with::

    conda install jupyter

See the `Jupyter documentation`_ for complete information about Jupyter notebooks.

__ Jupyter_

Fluxpart documentation and examples are here:

    `tutorial <https://fluxpart.readthedocs.io/en/latest/tutorial.html>`_

    `notebooks <https://github.com/usda-ars-ussl/fluxpart/tree/master/docs/notebooks>`_

Some example data files are here:

    `data <https://github.com/usda-ars-ussl/fluxpart/tree/master/tests/data>`_

If at some point it is desired to deactivate the **fluxpart**  environment,
then::

    conda deactivate
