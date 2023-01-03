Brief Introduction to Inputs
======================================
The dictionary ``CanteraTorchProperties`` is the original dictionary of DeepFlame. It reads in network related parameters and configurations. It typically looks like:

.. code-block::

    chemistry           on;
    CanteraMechanismFile "ES80_H2-7-16.yaml";
    transportModel "Mix";
    odeCoeffs
    {
        "relTol"   1e-15;
        "absTol"   1e-24;
    }
    inertSpecie        "N2";
    zeroDReactor
    {
        constantProperty "pressure";
    }

    splittingStretagy false;

    TorchSettings
    {
        torch on;
        GPU   off;
        log  on;
        torchModel "HE04_Hydrogen_ESH2_GMS_sub_20221101"; 
        coresPerNode 4;

    }
    loadbalancing
    {
            active  false;
            //log   true;
    }

In the above example, the meanings of the parameters are:

* ``CanteraMechanismFile``: the name of the reaction mechanism file.
* ``transportModel``: the default model is *Mix*, but other models including *UnityLewis* and *Multi* are also availabile.
* ``constantProperty``: property set to be constant during reaction. It can be set to *pressure* or *volume*.
* ``odeCoeffs``: the ode tolerance. 1e-15 and 1e-24 are used for network training, so they should be kept the same when comparing results with and without DNN. Default values are 1e-9 and 1e-15.
* ``TorchSettings``: all paramenters regarding the usage of DNN. This section will not be read in CVODE cases.
* ``torch``: the switch used to control the on and off of DNN. If users are running CVODE, this needs to be switched off.
* ``GPU``: the switch used to control whether GPU or CPU is used to carry out inference.
* ``torchModel``: name of network.     
* ``coresPerNode``: If you are using one node on a cluster or using your own PC, set this parameter to the actual number of cores used to run the task. If you are using more than one node on a cluster, set this parameter the total number of cores on one node. The number of GPUs used is auto-detected.



