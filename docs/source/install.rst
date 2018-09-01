=============
 Installation
=============

.. _Python: https://www.python.org
.. _NumPy: http://www.numpy.org
.. _SciPy: http://www.scipy.org
.. _Matplotlib: http://matplotlib.org
.. _IPython: https://ipython.org
.. _Jupyter: http://jupyter.org
.. _Jupyter documentation: http://jupyter.readthedocs.io/en/latest/index.html
.. _Pandas: https://pandas.pydata.org
.. _PyWavelets: http://pywavelets.readthedocs.io/en/latest
.. _miniconda: http://conda.io/miniconda.html
.. _anaconda: https://www.anaconda.com/download
.. _conda: http://conda.pydata.org/docs/user-guide/index.html
.. _guide: https://conda.io/docs/user-guide/getting-started.html

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

Recommended
~~~~~~~~~~~

* IPython_
* Jupyter_

Install
~~~~~~~

Recommended:

    Linux, Mac OS X::

        conda install anaconda-client
        conda env create ussl/fluxpart-env
        source activate fluxpart-env

    Windows::

        conda install anaconda-client
        conda env create ussl/fluxpart-env
        activate fluxpart-env

Alternatively:

    ::

        conda install -c ussl fluxpart

    or:

    ::

        pip install fluxpart


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
"Anaconda Promt". At the command prompt, enter::

    conda --version

Conda should respond with the version number of the installed package (e.g.,
``conda 4.5.5``). If you have an older, previously installed version of conda,
it is advised to update it  with::

    conda update conda

Next, install ``anaconda-client`` into the root environment::

    conda install anaconda-client

Finally, create a conda environment containing **fluxpart** and its
dependencies, as well as some useful libraries for plotting and analysis
(Jupyter_), by entering::

    conda env create ussl/fp020

To activate the **fluxpart** environment, enter:

    Linux, Mac OS X:::

        source activate fp020

    Windows::

        activate fp020

The command line prompt should now be prepended with ``(fp020)``,
indicating that the **fluxpart** environment is active in the shell session.

If you are new(ish) to Python_, a good tool for learning and interactively
building-up **fluxpart** analyses is the `Jupyter notebook`__. From the command
line, make and cd into a new working directory, e.g.::

    mkdir fluxnb
    cd fluxnb

An example notebook and high-frequency eddy covariance data file can be
downloaded with::

    anaconda download ussl/tutorial-notebook
    anaconda download ussl/tutorial-data

The Jupyter notebook application can be launched from the command line of an
active **fluxpart** shell session::

    jupyter notebook

The Jupyter dashboard will start in a web browser window, and look something
like this:

TODO

..
   .. image:: screenshot_jupyter_dashboard.png

Clicking on the ``tutorial-notebook.ipynb`` link will open the notebook in a
new browser tab:

TODO

..
   .. image:: screenshot_jupyter_notebook.png

Selecting a
cell in the notebook and hitting <Shift><Enter> executes the code in the cell.
See the `Jupyter documentation`_ for
complete information about Jupyter notebooks, and :ref:`fluxpart-tutorial` for
getting started with **fluxpart**.

__ Jupyter_

If at some point it is desired to deactivate the **fluxpart**  environment,
then:

    Linux, Mac OS X::

        source deactivate

    Windows::

        deactivate
