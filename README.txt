The Age-Layered Population Structure (ALPS) paradigm is a novel
metaheuristic for overcoming premature convergence by running multiple
instances of a search algorithm simultaneously.

In the field of optimization algorithms, one type of algorithmic
improvement is that of increasing the speed at which problems of
solvable difficulty can be solved. This is shown in papers in which
the comparison is on which algorithm can find the global optima of a
benchmark problem in the fewest number of evaluations. Another type of
algorithmic improvement is increasing the robustness of the
algorithm. That is, either increasing the reliability of finding the
global optima or being able to find better results then other
algorithms.

A typical approach to addressing premature convergence is to restart
the search algorithm, but one challenge is in deciding when the
population is truly stuck. Another problem is that all the information
learned from one run is not passed to the next run. An alternative to
restarting the entire search algorithm is to run multiple EAs
simultaneously and only restart one of them, and this is what is done
with ALPS.

The Age-Layered Population Structure (ALPS) paradigm is a novel
metaheuristic for overcoming premature convergence by running multiple
instances of a search algorithm simultaneously. A novel measure of age
is used to segregate individuals into different age-layers and then,
at regular intervals, the youngest layer is replaced with randomly
generated individuals.

ALPS was developed to be a way to make other search algorithms more
robust, especially on hard problems, but at the cost of potentially
slowing them down on easier problems.  For more details refer to the
web page or publications:
 http://idesign.ucsc.edu/projects/alps.html

Hornby, G. S. (2009) "Steady-State ALPS for Real-Valued Problems",
Proc. of the Genetic and Evolutionary Computation Conference, ACM
Press.

Hornby, G. S. (2009) "A Steady-State Version of the Age-Layered
Population Structure EA" , Genetic Programming Theory & Practice VII.

Hornby, G. S. (2006) "ALPS: The Age-Layered Population Structure for
Reducing the Problem of Premature Convergence", Proc. of the Genetic
and Evolutionary Computation Conference, ACM Press.

