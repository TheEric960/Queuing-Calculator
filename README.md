# Queuing-Calculator
The formulas in this Python script implement the simple queueing models described
in Chapter 6 of Banks, Carson, Nelson and Nicol, *Discrete-Event System 
Simulation*, 5th edition.

Supported Queues:
* M/G/1
* M/M/c
* M/G/c
* M/M/c/N
* M/M/c/K/K

Input Parameter Definitions:
 * `lmda`: "λ" i.e. arrival rate
 * `mu`: "μ" i.e. service rate
 * `c`: number of servers
 * `sigma2`: "σ<sup>2</sup>" i.e. variance of service time
 * `n`: system capacity, including customers in service
 * `k`: size of calling population

Output Parameter Definitions:
* `rho`: "ρ" i.e. utilization
* `l`: mean number in system
* `w`: mean time in system
* `wq`: mean time in queue
* `lq`: mean number in queue
* `p0`: probability of empty system
* `lmda_effective`: "λ<sub>effective</sub>" i.e. effective arrival rate

Usage Examples:

After `import queueing as q`:
 * `rho, l, w, wq, lq, p0 = q.eval_MG1(lmda=1.125, mu=2.35, sigma2=0.2)`
 * `rho, l, w, wq, lq, p0 = q.eval_MMc(lmda=3.6, mu=2.15, c=3)`
 * `rho, l, w, wq, lq = q.eval_MGc(lmda=5.42, mu=2.18, c=3, sigma2=0.56)`
 * `rho, l, w, wq, lq, p0, pN, lmda_effective = q.eval_MMcN(lmda=12.98, mu=3.47, c=4, n=15)`
 * `rho, l, w, wq, lq, p0, lmda_effective = q.eval_MMcK(lmda=2.65, mu=1.2, c=5, k=6)`

