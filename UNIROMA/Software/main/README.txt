3.2	Code
The code was developed with the following main scripts:
•	main.py
•	test_main.py
•	WX_Y.py
•	threshold_flexibility.py
Many functionalities are included in main.py which contains the following functions that are executed in series:
•	environment preparation
•	data gathering
•	state estimation
•	load flow
•	flexibility decision
•	calculation of flexibility effects
The first functionality environment preparation is developed in order to prepare the software to be executed properly. Therefore, all the necessary libraries are imported in first part of main.py. These libraries are necessary to execute the functions that will be described below. In particular, the script requires the following libraries:
•	py_dss_interface, this is necessary to interface the main script with the power flow calculator that is OpenDSS which is an open source software developed for studying the distribution networks.
•	Pygad, this library is necessary to implement the genetic algorithm within the optimization process (i.e., the state estimation and the calculation of flexibility thresholds).
•	Numpy, this library is necessary to make matrix calculations.
In order to properly execute the Python scripts it is necessary to set the environment in the related Python idle of the user of this software. In particular, it is necessary to set the right versions of these libraries. For the sake of completeness the right versions of these libraries are included in the online repository in which all the scripts are also available. In addition, in order to properly execute the py-dss-interface library, the library pandas is necessary. The related version of pandas library is also provided.
After the preparation of the environment, the additional functionality is the data gathering. This functionality was developed assuming that the microgrid or the assets are monitored by distributed sensors that can be managed by remote. In particular, it was decided to adopt MQTT broker as the interface between sensors and the broker. The functionality is implemented by means of a sequence of scripts which can be generally named WX_Y.py which include the following information:
•	W is the type of sensor
•	X is the number of the sensor
•	Y is the information that we want to get from the broker (e.g., active power or reactive power)
It is worth noting that these information are collected by distributed sensors that are normally installed in order to monitor some relevant nodes of the network (e.g., LV bus bar of a secondary substation, big power plant, big facility, line start). The implementation in the software is carried out so that these scripts are executed in parallel and therefore the information are received as synchronized as possible. In any case, delays among the data received that are lower than 5 seconds do not cause any problems to the developed services and the provided outputs are equally relevant. Therefore, delays among sensors are not considered and these errors are neglected assuming that all data are received simultaneously.
After the data gathering, the state estimation process start. For the sake of clarity the state estimation process is encapsulated in the main file of the software so that the process is transparent to the user and the developer that wants to modify the code. In general, the state estimation assumes that some nodes of the network are monitored while others are not. Given this circumstance, the state estimation calculates the missing values to be assigned to the nodes. The estimated values are provided according to some historical data and they are obtained through an optimization process that aims at minimizing the errors between the measured values and those calculated when the estimated values are applied. In this software, this optimization process is solved by a genetic algorithm. The usage of genetic algorithms allow to solve even those situations in which the number of monitored is lower that the minimum number of monitored nodes required for applying the traditional methods for state estimation.
The load flow calculation follows the state estimation. Indeed, the results of the state estimation are exploited for calculating the actual load flow in the network. After this step, the calculated values are definitely assigned to the nodes which are not monitored. The load flow calculation is exploited to assess some relevant Key Performance Indicators that are needed to activate the flexibility resources. In this current development, the indicator under consideration is the reverse power flow at the primary substation. According to the calculated reverse power flow, the thresholds for flexibility activation are set in a deterministic way as it is explained in the previous section.
After the calculations of the status of the network, the deterministic approach is exploited to identify the most suitable thresholds according to the status of the reverse power flow from the feeder to the transmission network. These values can be set or automatically calculated at the initialization process by means of threshold_flexibility.py. According to a certain level of reverse power flow, the identified values are stored for the next steps of the software execution.
The last step of the main script is the calculation of the effects of the flexibility on the network. Indeed, the results of the power flows and state estimation are modified considering the variations that would be caused by the flexibility sources. In particular, this part of the main script calculates again the key performance indicators in order to understand the effects of the flexibility resources in a distribution network. However the software is capable of considering the effects of the storage even assuming that the orders are actually enforced in the network (e.g., it can be assumed that some actuators are actually integrated with the main script). In this case, the comparison would be necessary to assess the effects on the network if the flexibility were not enforced.
As it is anticipated at the beginning of this section, some auxiliaries script were developed in order to aid the whole process. In particular, test_main.py was developed in order to parallelize the processes of data gathering. As it is explained such execution ensures that data are collected from the distributed sensors in a small timeframe so that the information delay do not cause any issues to the calculation (i.e., state estimation, flexibility enforcement). Another auxiliary script is WX_Y.py that is used to get the data from the field. The developed script are invoked in test_main.py.
The last auxiliary script to be described is threshold_flexibility.py. This script is useful to determine the thresholds for flexibility orders according to a specific status of the network (e.g., high reverse power flow). It implements the model of flexibility calculation previously described and it provides as output the thresholds for the flexibility enforcement.



Use case
The network section selected for the case study well represents a section is an MV feeder that supplies 5 secondary substations, as shown in the below figure. Therefore, the opendss model of the network is implemented in the main directory of the repository. The dss model consists of the following scripts:
•	Run.dss
•	Master.dss
•	Loads.dss
•	Lines.dss
•	Storage.dss
•	Transformers.dss
With respect to the data, some datasets are preloaded in the example whilst the connection to the broker is still not available since it is currently operated on a private network.
These are the functionalities implemented by the scripts: 
•	Run.dss executes the master in the proper directory 
•	Master.dss compiles all the files about loads, lines and other DN elements. Then, it provides the needed outuputs 
•	Loads.dss describes the load absorption, type of connection and the nodes to be connected in the network, 
•	Lines.dss includes the description of the lines from the topological perspective, as well as it contains information about the line parameters (i.e., resistance and reactance of the lines); additional information can be also included such as the line ampacity, 
•	Storage.dss contains information about the energy storage systems located in the DN. This script allows to take into account the state of charge of the batteries as well as the maximum power that can be absorbed or injected by the storages. 
•	Transformers.dss is similar to Lines.dss but it exploits the transformer object and it contains information about transformer parameters. 
