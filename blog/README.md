# reinforced genetic algorithm for structure-based drug design


This blog is written to introduce our recent NeurIPS 2022 paper [1]: Reinforced Genetic Algorithm for Structure-based Drug Design. 


The paper proposes a reinforcement learning-based genetic algorithm for structure-based drug design. 
The following figure illustrate the difference between genetic algorithm and reinforced genetic algorithm. 


<p align="center"><img src="fig/RGA.png" alt="logo" width="860px" /></p>



## Problem: structure-based Drug Design


Rapid drug discovery that requires less time and cost is of significant interest in pharmaceutical science. Structurebased drug design (SBDD) that leverages the three-dimensional (3D) structures of the diseaserelated proteins to design drug candidates is one primary approach to accelerate the drug discovery processes with
physical simulation and data-driven modeling. According to the lock and key model, the molecules that bind tighter to a disease target are more likely to expose bioactivity against the disease, which has been verified experimentally. 
As AlphaFold2 has provided accurate predictions to most human proteins [3], SBDD has a tremendous opportunity to discover new drugs for new targets that we cannot model before [2].



## genetic algorithm

Traditional combinatorial optimization methods such as genetic algorithms (GA) have demonstrated state-of-the-art performance in various molecular optimization tasks. However, they do not utilize protein target structure to inform design steps but rely on a random-walk-like exploration, which leads to unstable performance and no knowledge transfer between different tasks despite the similar binding physics.



## reinforced genetic algorithm


To achieve a more stable and efficient SBDD, we propose \fullname (\mname) that uses neural models to prioritize the profitable design steps and suppress random-walk behavior. 
The neural models take the 3D structure of the targets and ligands as inputs and are pre-trained using native complex structures to utilize the knowledge of the shared binding physics from different targets and then fine-tuned during optimization. 





## Experiment


**Optimization ability** 
We report the optimization performance of all the methods as follows. 
We evaluate all the methods on all targets and report each metric's mean and standard deviations across all targets. 
Arrows indicate the direction of better performance. 
For each metric, the best method is underlined and the top-3 methods are bolded. 
RGA-pretrain and RGA-KT are two variants of RGA that without pretraining and without training on different target proteins, respectively. 
Our result shows reinforced genetic algorithm (RGA) achieves the best performance in TOP-100/10/1 scores among all methods we compared. Compared to traditional genetic algorithm (Autogrow 4.0 [2]), RGA's better performance in docking score demonstrates that the policy networks contribute positively to the chemical space navigation and eventually help discover more potent binding molecules.


<p align="center"><img src="fig/results.png" alt="logo" width="740px" /></p>


**Sample efficiency** Especially in SBDD, when the oracle functions are expensive molecular simulations, robustness to random seeds is essential for improving the worst-case performance of algorithms. One of the major issues in traditional GAs is that they have a significant variance between multiple independent runs as they randomly select parents for crossover and mutation types. To examine this behavior, we run five independent runs for RGA, Autogrow 4.0 and graph-GA (three best baselines, all are GA methods) on all targets and plot the standard deviations between runs. As shown in the following figure, with policy networks guiding the action steps, we observed that the random-walk behavior in Autogrow 4.0 was suppressed in RGA, indicated by the smaller variance. Especially in the later learning phase (after 500 oracle calls), the policy networks are fine-tuned and guide the search more intelligently. This advantage leads to improved worst-case performance and a higher probability of successfully identifying bioactive drug candidates with constrained resources.



<p align="center"><img src="fig/sns.png" alt="logo" width="350px" /></p>



We have open-sourced our code in https://github.com/futianfan/reinforced-genetic-algorithm.

## References

[1] Tianfan Fu*, Wenhao Gao*, Connor W. Coley, Jimeng Sun. Reinforced Genetic Algorithm for Structure-based Drug Design. Neural Information Processing Systems (NeurIPS) 2022. 

[2] Spiegel, J. O. and Durrant, J. D. Autogrow4: an opensource genetic algorithm for de novo drug design and
lead optimization. Journal of cheminformatics, 12(1): 1–16, 2020.

[3] Jumper, John, et al. "Highly accurate protein structure prediction with AlphaFold." Nature 596.7873 (2021): 583-589. 


