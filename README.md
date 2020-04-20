# LeonardYM

LeonardYM is a software developed in C++ for Monte-Carlo simulations of four-dimensional Yang-Mills theories, ranging from QCD to supersymmetric models. The code is designed to allow many different theories to be simulated within the same framework, rather than being focused and optimized for a single specific target theory. Parameters such as the number of colors of the gauge group SU(N) or the representation of fermion fields can be freely chosen at compile time, allowing a large flexibility on the structure of the theories that can be simulated.

## Getting started

Details and instructions on how to run and compile the code can be found in the folder [documentation/](documentation/documentation.pdf). Example configuration scripts are available in the folder [configuration_scripts/](configuration_scripts/).

### Prerequisites

* Gauge-field variables are defined from the matrix library [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* High-level parallelism is based on OpenMP and MPI. 
* The code uses the library program-options of [boost](https://www.boost.org/)  

## Authors

* **Stefano Piemonte** - *Initial work* - [spiemonte](https://github.com/spiemonte)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the GNU GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* The development of the code has started during the work for the PhD thesis [N=1 supersymmetric Yang-Mills theory on the lattice](https://www.uni-muenster.de/imperia/md/content/physik_tp/theses/muenster/piemonte_dr.pdf). The software has been used mainly to perform Monte-Carlo simulations of supersymmetric Yang-Mills theories by the DESY-Muenster collaboration, whose support and help during the development of the code has been crucial. 

