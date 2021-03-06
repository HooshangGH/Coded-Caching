# Asynchronous Coded Caching
Here, I want to present overall idea of our work on asynchronous coded caching problem and explain the code we used to generate the simulation results in our work. For more detail about our work please take a look at <a href="https://arxiv.org/pdf/1907.06801.pdf" title="this">this</a> and <a href="https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8006967" title="this">this</a> papers.
## Original Coded Caching Problem
In their pioneering work, Maddah-Ali and Niesen considered the usage of coding in the caching problem. In this so-called "coded caching" setting, there is a server containing $N$ files, $K$ users each with a cache that can store up to $M$ files. The users are connected to the server via an error-free shared link.

The system operates in two distinct phases. 
* In the *placement phase* the content of the caches is populated by server. This phase does not depend on the future requests of the users which are assumed to be arbitrary. 
* In the *delivery phase* each user makes a request and the server transmits potentially coded signals to satisfy the requests of the users.

The work of Maddah-Ali and Niesen demonstrated that significant reductions in the network traffic were possible as compared to conventional caching. Crucially, these gains continue to hold even if the popularity of the files is not taken into account.


## Asynchronous Coded Caching
While this is a significant result, the original formulation of the coded caching problem assumes that the user requests are synchronized, i.e., all file requests from the users arrive at the server at the same time. Henceforth, we refer to this as the synchronous setting. From a practical perspective, it is important to consider the asynchronous setting where user requests arrive at different times.
	In this case, a simple strategy would be to wait for the last request to arrive and then apply the scheme of Maddah-Ali and Niesen. Such a strategy will be quite good in terms of the overall rate of transmission from the server. However, this may be quite bad for an end user's experience, e.g., the delay experienced by the users will essentially be dominated by the arrival time of the last request.
	
<p align="center">
  <img src="Netflix.jpg" width="300" height="250">
</p>

In this work we formulate and study the coded caching problem when the user requests arrive at different times. Each user has a specific deadline by which his/her demand needs to be satisfied. The goal is to schedule transmission of packets so that each user is able to recover the requested file from the transmitted packets and his/her cache content within the prescribed deadline. We present algorithms for both the offline and online versions of this problem.

### Problem Formulation
We assume that time $\tau \geq 0$ is slotted. Let $[n]$ denote the set $\{1, \ldots, n\}$ and the symbol $\oplus$ represent the XOR operation. We assume that the server contains $N\geq K$ files\footnote{We assume that $N\geq K$ as it corresponds to the worst case rate where each of the $K$ users can request a different file. Furthermore, it is also the more practical scenario.} denoted by $W_{n}, n = 1, \dots, N$. The subfiles are denoted by $W_{n,f}$ so that $W_n = \{W_{n,f}: f \in [F]\}$ and the cache of user $i$ by $Z_i \subseteq \{ W_{n,f}: \ n \in [N], \ f \in [F] \}$. $Z_i$ contains at most $M_iF$ subfiles. In the delivery phase, user $i$ requests file $W_{d_i}$, where $d_i \in [N]$, from the server.
We let $\Omega^{(i)}$ denote the indices of the subfiles that are not present in the $i$-th user's cache, i.e.,
$$
\Omega^{(i)} = \{ f:\ f \in [F], \ W_{d_i,f} \notin Z_i \}.
$$
The equations in the delivery phase are assumed to be of the *all-but-one* type.
*All-but-one equation*. Consider an equation $E$ such that
$$
E = \oplus_{l=1}^\ell W_{d_{i_l}, f_{l}}.
$$
We say that $E$ is of the all-but-one type if for each $l \in [\ell]$, we have $W_{d_{i_l}, f_{l}} \notin Z_{i_l}$ and $W_{d_{i_l}, f_{l}} \in Z_{i_{k}}$ for all $k \in [\ell] \setminus \{l\}$.
\end{definition}
It is evident that an all-but-one equation transmitted from the server allows each of the users participating in the equation to recover a missing subfile that they need. The asynchronous coded caching problem can be formulated as follows.
	
*Inputs.*
* *User requests.* User $i$ requests file $W_{d_i}$, with $d_i \in [N]$ at time $T_i$. %User $i$'s request arrives at the server at time $T_i$.
* *Deadlines.* The $i$-th user needs to be satisfied by time $T_i + \Delta_i$, where $\Delta_i$ is a positive integer.
* *Transmission delay.* Each subfile needs $r$ time-slots to be transmitted over the shared link, i.e., each subfile can be treated as equivalent to $r$ packets, where each packet can be transmitted in one time slot.

As the problem is symmetric with respect to users, w.l.o.g. we assume that $T_1 \leq T_2 \leq \ldots \leq T_K$. Let $T_{\max} = \max_i (T_i + \Delta_i)$. Note that upon sorting the set of arrival times and deadlines, i.e., $\cup_{i=1}^{K} \{T_i, T_i + \Delta_i\}$, we can divide the interval $[T_1, T_{\max})$ into *at most* $2K-1$ non-overlapping intervals. Let the integer $\beta$, where $1 \leq \beta \leq 2K-1$ denote the number of intervals.
	

Thus, $U_\ell$ is the set of active users in time interval $\Pi_\ell$ and $D_\ell$ is the corresponding set of active file requests.
	
*Outputs.*
* *Transmissions at each time slot.* If the problem is feasible, the schedule specifies which equations (of the all-but-one type) need to be transmitted at each time. The schedule is such that each user can recover all its missing subfiles within its deadline. The equations transmitted at time $\tau \in \Pi_\ell$ only depend on $D_\ell$. 

### Offline Case
In this section, we discuss the offline version of the problem where the server has the knowledge of the arrival times/deadlines of all the requests at $\tau = 0$. We have written a python class (OffACC.py) to return offline solution for a problem instance. You can find the code <a href="https://github.com/HooshangGH/Coded-Caching/">here</a> and download it from <a href= "OffACC.py">here</a>. For simplicity, each subfile $W_{n,f}$ is denoted by index $j = (n-1)*F + f -1$. The class must be initialized by the following parameters:
```js
def __init__(self, N, K, M, F, arrival_times, deadlines, cache_contents, Omega):
        '''
        :param N:               number of files
        :param K:               number of users
        :param M:               cache size normalized by file size
        :param F:               size of each file
        :param arrival_times:   arrival times of the users
        :param deadlines:       deadline of each user
        :param cache_contents:  cache content of each user
	'''
```
Basically 
* arrival_times $= [T_1,\ldots, T_K]$ and deadlines $$=[T_1+\Delta_1,\ldots, T_K+\Delta+K]$$. 
* Also, cache_contents is an array of array of size $K$ so that cache_contents[i] = $$[(n-1)*N+f-1 ~for ~ W_{n,f} \in Z_i]$$. * Similarly, Omega is an array of array so that Omega[i] = $$[(d_i-1)*N+f-1 ~for~ W_{d_i,f} \in \Omega^{(i)}]$$.

To solve the offline problem, you need to call ```make_n_solve_LP``` function after creating an object of the class and initializing it. Output to this function is solution of the LP and rate (minimum number of transmitted equations (packets)). If the problem is feasible, the rate will be -1. Here is an example for Example 1 in <a href="https://arxiv.org/pdf/1907.06801.pdf">our paper</a>. 
```js
from OffACC import *
# setting of the problem
N = 3
K = 3
M = 1
F = 3
arrival_times = [0, 1, 3]
deadlines = [5, 5, 5]
cache_contents = [[3, 4, 8], [0, 5, 6], [3, 4, 7]]
Omega = [[0, 1, 2], [3, 4], [6, 8]]
# object of the class
objOff = OffACC(N, K, M, F, arrival_times, deadlines, cache_contents, Omega)
# calling make_nsolve_LP
off_solution, rate = objOff.make_n_solve_LP()
print(off_solution)
print(rate)
```
Also, here is the output when you run this code
```js
{'{0}-1': 1.0, '{0}-0': 1.0, '{2}-2': 1.0, '{1,2}-2': 1.0, '{0,1}-1': 1.0}
5.0
```
```off_solution``` is set of nonzero $x_U(\ell)$'s after solving LP in (1) in <a href="https://arxiv.org/pdf/1907.06801.pdf">our paper</a>.

### Online Case
In the online scenario, at time $\tau$ only information about the already arrived requests are known to the server, i.e., it only knows $T_i$, $d_i$ and $\Delta_i$ for $i \in [K]$ such that $T_i \leq \tau$.

We have ```Decentralized_CC``` in our ```ODCC.py``` file that can be used for the online case. You can download it from <a href="ODCC.py">here</a>. The class has to be initialized in this way
```js
def __init__(self, K, N, F, n, r,
                 cache_contents,
                 requested_files,
                 arrival_times,
                 deadlines,
                 threshold_eta):
		 
```
The only differences from the offline case are that we have $n$, $r$, ```requested_files```, and ```threshold_eta```. We have talked about these parameters in <a href="https://arxiv.org/pdf/1907.06801.pdf">our paper</a>. Unlike offline case, instead of feeding $\Omega^{(i)}$'s to the class, we use ```requested_files``` $=[d_1, \ldots, d_K]$. Here we show how we use our code for the online case of Example 1 in <a href="https://arxiv.org/pdf/1907.06801.pdf">our paper</a>.

```js
from OffACC import *
from ODCC import *
# setting of the problem
N = 3
K = 3
M = 1
F = 3
arrival_times = [0, 1, 3]
deadlines = [5, 5, 5]
cache_contents = [[3, 4, 8], [0, 5, 6], [3, 4, 7]]

requested_files = [i+1 for i in range(K)]

eta = 0.0
n = 1
r = 1

objOnline = Decentralized_CC(K, N, F, n, r, cache_contents, requested_files, arrival_times, deadlines, eta)

rate, online_solution = objOnline.run_algo_mod_eta()

print(rate)
print(online_solution)
```
When running this code, the following will be printed out.
```
-1
None
```
This means that our algorithm was not able to come up with a solution that satisfies all requestes within their deadline. The reason is that time $\tau=2$ the algorithm prefers serving user 2 instead of the user 1 and since use 1 and user 2 can't simultanouesly benefit from a *all-but-one* type of euation thus the remaining time after $\tau=3$ is not enough to transmit 3 missing subfiles of users 1 and 2 only in 2 time slots.


## Dual Decomposition
The complexity of the solving the LP does grow quite quickly (cubic) in the problem parameters. The LP in (1) in <a href="https://arxiv.org/pdf/1907.06801.pdf">our paper</a> can however be modified slightly so that the corresponding dual function is such that it can be evaluated by solving a set of *decoupled minimum cost network flow optimizations*. Please take a look at this paper for more detail.

I also wrote a cpp code to solve LP in (1) in <a href="https://arxiv.org/pdf/1907.06801.pdf">our paper</a> by dual decomposition idea. The code is available <a hfre="MCFgraph">here</a>. You must creat an object of ```ProblemInstance``` class and run it. I used <a href="https://lemon.cs.elte.hu/trac/lemon/wiki/MinCostFlowData">Lemon Package</a> for solving minimum cost flow problems. You need to install Lemon package for using my code. Here is how I use this code.
```js
const int K = 20;
const int t = 2;
const int r = 1;
int Tmaxitr = 1000;
double alpha = 0.99;
int offItrs = 20;
// Arrival times and deadlines
vector<int> deadline = {0, 10, 20, 79, 83, 108};
vector<int> arrival = {51, 126, 91, 138, 162, 220};
// creating an object of the class
ProblemInstance prb_instc(arrival, deadline, K, t, r);
// Constructing the problem instance and associated MCF graph
prb_instc.ConstructProblem();
// Solving the problem instance.
prb_instc.SolveProblem(Tmaxitr, alpha, offItrs);

```

Please feel free to reach out to me about the codes.
