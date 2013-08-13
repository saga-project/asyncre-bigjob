ASyncRE-BigJob
==============

ASynchronous Replica Exchange (ASyncRE-BigJob) is an extensible Python package enabling file-based asynchronous parallel replica exchange molecular simulations using the BigJob distributed computing infrastructure.

Replica Exchange (RE) is a popular generalized ensemble approach for the efficient sampling of conformations of molecular systems. In RE, the system is simulated at several states differing in thermodynamic environmental parameters (temperature, for example) and/or potential energy settings (biasing potentials, etc). Multiple copies (replicas) of the system are simulated at each state in such a way that, in addition to traveling in conformational space, they also travel in state space by means of periodic reassignments of states to replicas. Traditional synchronous implementations of RE are limited in terms of robustness and scaling because all of the replicas are simulated at the same time and state reassignments require stopping all of the replicas. In Asynchronous RE replicas run independently from each other, allowing simulations involving hundreds of replicas on distributed, dynamic and/or unreliable computing resources.

In ASyncRE, the BigJob framework is used for launching, monitoring, and managing replicas. State exchanges are performed for idle replicas via the filesystem by extracting and modifying data on the input/output files of the MD engine while other replicas continue to run. Support for arbitrary RE approaches is provided by simple user-provided adaptor modules which in general do not require source code-level modifications of legacy simulation engines. Currently, adaptor modules exist for umbrella sampling RE with the AMBER engine and BEDAM binding free energy calculations with IMPACT.

Web Pages & Mailing Lists
-------------------------

ASyncRE-BigJob: https://github.com/saga-project/asyncre-bigjob/wiki

BigJob: https://github.com/saga-project/BigJob/wiki

async-replica-exchange@googlegroups.com

async-replica-exchange-devel@googlegroups.com

Installation
------------

The ASyncRE packages primarily depends on BigJob. To install it follow the instructions on the BigJob site:  http://http://saga-project.github.io/BigJob/. Note that BigJob requires a working queuing system and a redis server. ASyncRE also depends numpy and configobj, which are easily installed from PiP: 

    pip install numpy
    pip install configobj

ASyncRE is currently distributed only by git:

    git clone https://github.com/saga-project/asyncre-bigjob.git
    cd asyncre-bigjob
    python setup.py install

A distribution archive can be created by issuing the command:

    python setup.py sdist

after which async_re-<version>.tar.gz will be found under dist/

Installation from the distribution archive:

    cd dist
    pip install async_re-<version>.tar.gz


Test
----

To test execute the "date" application

    python date_async_re.py command.inp

which will spawn a bunch of /bin/date replicas.

See additional sample application files under the examples/ subdirectory.

Documentation
-------------

Documentation is in the doc/ subdirectory and on the ASyncRE wiki (https://github.com/saga-project/asyncre-bigjob/wiki, currently under construction).

