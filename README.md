dcsxx-cloud-gt
==============

Game-theoretic approach to coalition formation in green cloud federations.

This code has been used in the following paper:

> Marco Guazzone, Cosimo Anglano and Matteo Sereno.
>
> *A Game-Theoretic Approach to Coalition Formation in Green Cloud Federations*
>
> Proc. of the 14th IEEE/ACM International Symposium on Cluster, Cloud and Grid Computing (CCGrid), pp. 618-625, 2014.
>
> doi:[10.1109/CCGrid.2014.37](http://dx.doi.org/10.1109/CCGrid.2014.37).

Please, cite this code as follow (BibTex format):

	@INPROCEEDINGS{Guazzone-2014-Game, 
		author = {Marco Guazzone and Cosimo Anglano and Matteo Sereno}, 
		booktitle = {Cluster, Cloud and Grid Computing (CCGrid), 2014 14th IEEE/ACM International Symposium on}, 
		title = {A Game-Theoretic Approach to Coalition Formation in Green Cloud Federations}, 
		year = {2014}, 
		pages = {618--625}, 
		keywords = {cloud federation; cooperative game theory; coalition formation; cloud computing; green computing}, 
		doi = {10.1109/CCGrid.2014.37}, 
		month = {May},
	}

Overview
--------

Federations among sets of Cloud Providers (CPs), whereby a set of CPs agree to mutually use their own resources to run the Virtual Machines (VMs) of other CPs, are considered a promising solution to the problem of reducing the energy cost.

In (Guazzone,2014), we address the problem of federation formation for a set of CPs, whose solution is necessary to exploit the potential of cloud federations for the reduction of the energy bill.

 We devise an algorithm, based on cooperative game theory, that can be readily implemented in a distributed fashion, and that allows a set of CPs to cooperatively set up their federations in such a way that their individual profit is increased with respect to the case in which they work in isolation.

We show that, by using our algorithm and the proposed CPs' utility function, they are able to self-organize into Nash-stable federations and, by means of iterated executions, to adapt themselves to environmental changes.

Numerical results are presented to demonstrate the effectiveness of the proposed algorithm.


Requirements
------------

* [Boost](http://www.boost.org) >= v1.55
* [dcsxx-commons](https://github.org/sguazt/dcsxx-commons) >= v2
* [gtpack](https://github.com/sguazt/gtpack) v1
* [IBM ILOG CPLEX Optimizer](http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/index.html) >= v12.6


Compiling
----------

1. Edit the provided `Makefile` to change configuration parameters
2. Run `make clean all`


Running
-------

Usage:
```
./src/sim <options>
```

The program accepts the following command line options:

- `--csv <file>`
	<br/>
	Export all the analyzed coalition onto a CSV file `<file>`.
- `--formation {'merge-split'|'nash'|'pareto'|'social'}`
	<br/>
	The coalition formation strategy. Can be one of the following:
	- `'merge-split'`: to form Merge/split-stable partitions
	- `'nash'`: to form Nash-stable partitions
	- `'pareto'`: to form Pareto-optimal partitions
	- `'social'`: to form social-optimum partitions
- `--help`
	<br/>
	Show a help message.
- `--opt-relgap <num>`
	<br/>
	A real number `<num>` in [0,1] used to set the relative gap parameter of the optimal solver.
- `--opt-tilim <num>`
	<br/>
	A real positive number `<num>` used to set the maximum number of seconds to wait for the termination of the optimal solver.
- `--payoff {'banzhaf'|'norm-banzhaf'|'shapley'}`
	<br/>
	The coalition value division strategy. Can be one of the following:
	- `'banzhaf'`: the Banzhaf value
	- `'norm-banzhaf'`: the normalized Banzhaf value
	- `'shapley'`: the Shapley value
- `--rnd-genvms`
	<br/>
	Enable the random generation of the number of VMs for each CP (see sets <code>J<sub>i</sub></code> in (Guazzone,2014)).
	The number of VMs is randomly generated according to a Uniform(0,Nij) discrete probability distribution, where Nij is the maximum number of class-j VMs owned by CP i as defined by the `cip_num_vms` property of the scenario file.
	If provided, this option overrides the `cip_num_vms` property of the scenario file, which is now interpreted as the maximum number of VMs owned by each CP.
- `--rnd-genpmsonoff`
	<br/>
	Enable the random generation of host power states for each CP and host.
	The PM power state is randomly generated according to a Bernoulli(0.5) discrete probability distribution.
	If provided, this option overrides the `cip_pm_power_states` property of the scenario file.
- `--rnd-genpmsonoffcosts`
	<br/>
	Enable the random generation of host switch-on/off transition costs for each CP and host class (see parameters <code>L<sub>k</sub>*E_k</code> and <code>S<sub>k</sub>*E_k</code> in (Guazzone,2014)).
	The switch-on/off costs are computed as the product of the electricity price, the maximum host power consumption and the time taken to complete the switch-on or switch-off operation.
	Host switch-on/off time is randomly generated according to a Normal probability distribution, by following the characterization found in (Meisner,2009).
	If provided, this option overrides the `cip_pm_asleep_costs` and `cip_pm_awake_costs` properties of the scenario file.
- `--rnd-genvmsmigrcosts`
	<br/>
	Enable the random generation of CP-to-CP VM migration costs for each CP and VM class (see parameters <code>G<sub>c1,c2,k</sub></code> in (Guazzone,2014)).
	The migration cost is computed as the product between the data transfer cost rate and the data size to transfer, and assuming that our algorithm activates every 12 hours.
	The data transfer cost rate is taken from Amazon EC2 pricing (Amazon,2013).
	The data size to transfer is computed by assuming that data are persistently transferred during the migration time at a fixed data rate of 100 Mbit/sec for all CPs.
	The migration time is randomly generated according to a Normal distribution whose parameters are obtained from (Akoush,2010).
	If provided, this option overrides the `cip_to_cip_vm_migration_costs` property of the scenario file.
- `--rnd-numit <num>`
	<br/>
	Set the number `<num>` of times that the given scenario must be run.
	Each run generates new random values for the parameters controlled by the `--rnd-gen*` command-line options (if specified).
- `--rnd-seed <num>`
	<br/>
	Set the seed `<num>` to use for random number generation.
- `--scenario <file>`
	<br/>
	The path to the file `<file>` describing the scenario to use for the experiment.

You can find a sample scenario file in the `examples` directory.


References
----------
- **(Akoush,2010)** S. Akoush, R. Sohan, A. Rice, A.W. Moore and Andy Hopper.
*Predicting the Performance of Virtual Machine Migration*,
Proc. of the 2010 IEEE International Symposium on Modeling, Analysis and Simulation of Computer and Telecommunication Systems (MASCOTS), pp. 37-46, 2010.
- **(Amazon,2013)** Amazon EC2.
*Amazon EC2 Data Transfer Pricing*,
2013, Online: https://aws.amazon.com/ec2/pricing/#DataTransfer
- **(Guazzone,2014)** M. Guazzone, C. Anglano and M. Sereno.
*A Game-Theoretic Approach to Coalition Formation in Green Cloud Federations*,
Proc. of the 14th IEEE/ACM International Symposium on Cluster, Cloud and Grid Computing (CCGrid), pp. 618-625, 2014.
- **(Meisner,2009)** D. Meisner, B. Gold and T. Wenisch.
*PowerNap: Eliminating Server Idle Power*,
Proc. of the 14th International Conference on Architectural Support for Programming Languages and Operating Systems (ASPLOS), pp. 205-216, 2009.
